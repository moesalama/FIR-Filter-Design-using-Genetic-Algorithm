close all
clear
clc

% Define filter specifications
Fp = 0.4;  % Passband frequency (normalized, 0 to 1)
Fs = 0.6;  % Stopband frequency (normalized, 0 to 1)
Fc = (Fp+Fs)/2; % Cutoff frequency (normalized, 0 to 1)
Rp = 0.1;  % Passband ripple (dB)
As = 60;   % Stopband attenuation (dB)
wp = Fp*pi;
ws = Fs*pi;
wc = Fc*pi;
delta_p = (10^(Rp/20)-1)/(10^(Rp/20)+1);
delta_s = (1+delta_p)*10^(-As/20);
delta_w = ws-wp;
M = ceil((-10*log10(delta_p*delta_s)-13)/(2.285*delta_w));

% Adjust the order to ensure it's even
if mod(M, 2) ~= 0
    M = M + 1;
end

% Calculate ideal sinc coefficients for a low-pass filter
alpha = (M - 1) / 2;
n = 0:M-1;
h = Fc*sinc(Fc * (n - alpha));
ideal_coefficients = h;

%% Run this for 20 different seeds and record values for statistical inference
record = {};
for k = 1:20
    rng(k);
    % Fitness function handle
    fitnessFcn = @(coefficients) fir_fitness_function(coefficients, ideal_coefficients, Fp, Fs, Rp, As, M);

    % Genetic Algorithm options
    options = optimoptions('ga', ...
        'PopulationSize', 100, ...
        'MaxGenerations', 600, ...
        'SelectionFcn', @selectiontournament, ... % Tournament selection
        'CrossoverFcn', @crossoverarithmetic, ... % Arithmetic crossover
        'MutationFcn', @mutationadaptfeasible, ... % Adaptive feasible mutation
        'CrossoverFraction', 0.4, ... % Crossover fraction
        'Display', 'iter');

    % Run Genetic Algorithm
    [optimal_half_coefficients, optimal_fitness] = ga(fitnessFcn, (M/2)+1, [], [], [], [], [], [], [], options);

    % Construct the full set of coefficients for Type-1 FIR filter
    optimal_coefficients = [optimal_half_coefficients, fliplr(optimal_half_coefficients(1:end-2))];


    % Display results
    disp('Optimal Coefficients:');
    disp(optimal_coefficients);
    disp('Optimal Fitness:');
    disp(optimal_fitness);
    record{1,k} = optimal_coefficients;
    record{2,k} = optimal_fitness;
end

% Calculate the frequency response of the windowed filter
[H_optimal, f] = freqz(optimal_coefficients, 1, 1024, 'whole');
H_optimal_dB = 20*log10(abs(H_optimal)/max(abs(H_optimal)));
f = f / pi; % Normalize frequency

% Plot the frequency response
figure;
plot(f, H_optimal_dB);
title('Frequency Response of the Optimized Lowpass FIR Filter');
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Magnitude (dB)');
grid on;
hold on;

% Plot passband and stopband constraints
[H_ideal,f1] = freqz(ideal_coefficients,1,1024,"whole");
H_ideal_dB = 20*log10(abs(H_ideal)/max(abs(H_ideal)));
f1 = f1/pi;
plot(f1,H_ideal_dB)
legend('Frequency Response','Ideal Response');
%ylim([min(min(abs(H_ideal-dB)),min(abs(H_optimal_dB))) 1])
hold off;
err = h-optimal_coefficients;
figure
stem(err)

coeffs = [];
for k = 1:20
    coeffs = [coeffs;record{1,k}];
end

avgCoeff = [];
for k = 1:28
    avgCoeff(k) = mean(coeffs(:,k));
end

fitnessz = [];
for k = 1:20
    fitnessz = [fitnessz;record{2,k}];
end
