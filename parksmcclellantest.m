close all
clear
clc
%%
delta_p = 0.05;     % Passband error
delta_s = 0.01;     % Stopband error
wp = 0.4*pi;        % Passband frequency
ws = 0.6*pi;        % Stopband frequency
wc = (wp+ws)/2;     % CUtoff frequency
M = 1+ceil((-10*log10(delta_p*delta_s)-13)/(2.324*(ws-wp)));
% Filter order
w = [0, wp, ws, pi]; % Critical requencies
w_by_pi = w/pi;     % Normalized axis
A = [1 1 0 0];      % Critical amplitude factors
W = [1 delta_p/delta_s];    % Weighting function

h = firpm(M,w_by_pi,A,W);   % Parks-McClellan algorithm filter coefficients
n = 0:M;            % Discrete time vector
stem(n,h)           % Plot h[n]
xlabel('Time $n$','Interpreter','latex')
ylabel('Impulse response $h[n]$','Interpreter','latex')
[H,omega] = freqz(h,1,512);

% Ideal filter coefficients
Fc = wc/pi;
alpha = (M - 1) / 2;
n_ideal = 0:M;
h_ideal = Fc*sinc(Fc * (n_ideal - alpha));
[H_ideal,f] = freqz(h_ideal,1,512);

% Indices for critical frequencies in the omega vector
i1 = find(abs(omega - wp) <= 0.005);
i2 = find(abs(omega - ws) <= 0.005);
i3 = find(abs(omega-(wp+ws)/2)<=0.005);

% Expanded weighting function (size match)
W2 = [ones(i1(1),1); zeros(i2(1)-i1(1),1); ones(512-i2(1),1)];
Hd = [ones(i1(1),1); zeros(512-i1(1),1)];

% Error (objective) function to be minimized (for the algorithm)
E = abs(W2.*(abs(H)-abs(Hd)));
figure
plot(omega/pi,20*log10(abs(H)))
xlabel('Normalized Frequency $\omega$ ($\times \pi$ rad/sample)','Interpreter','latex')
ylabel('Gain (dB)','Interpreter','latex')
ylim([-80 5])
hold on
plot(f/pi,20*log10(abs(H_ideal)))
hold off

figure
plot(omega/pi,E,omega/pi,delta_p*ones(size(omega)),'k--',omega/pi,delta_s*ones(size(omega)),'k--')
ylim([0 0.06])
xlabel('Normalized Frequency $\omega$  ($\times\pi$ rad/sample)','Interpreter','latex')
ylabel('Error, $E(\omega)$','Interpreter','latex')