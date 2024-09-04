close all
clear
clc
%%
wp = 0.4*pi;        % Passband frequency
ws = 0.6*pi;        % Stopband frequency 
wc = 0.5*pi;        % Cutoff frequency
delta_w = ws-wp;    % Transition width
M = ceil(4*pi/delta_w); % Filter coefficients formula (rectangular window)
% Prepare practical filter
n = 0:M-1;
w = boxcar(M)';
hd = (wc/pi)*sinc((wc/pi)*(n-M/2));
h = hd.*w;
% Frequency response
[H,f] = freqz(h,1,512);
[Hd,f] = freqz(hd,1,512);
H_dB = 20*log10(abs(H)/max(abs(H)));
Hd_dB = 20*log10(abs(Hd)/max(abs(Hd)));
subplot(2,1,1)
plot(f/pi,H_dB)
xlabel('Normalized Frequency $\omega$ ($\times \pi$ rad/sample)','Interpreter','latex')
ylabel('Gain (dB)','Interpreter','latex')
hold on
plot(f/pi,-21*ones(size(f)),'--k')
plot(f/pi,Hd_dB,'r');
ylim([min(H_dB) 5])
grid on
subplot(2,1,2)
plot(f/pi,angle(H))
xlabel('Normalized Frequency $\omega$ ($\times \pi$ rad/sample)','Interpreter','latex')
ylabel('Phase (radians)','Interpreter','latex')
grid on