function fitness = fir_fitness_function(coefficients, ideal_coefficients, Fp, Fs, Rp, As, N)
    % Construct the full set of coefficients for Type-1 FIR filter
    % Adjust the order to ensure it's even
    if mod(N, 2) ~= 0
        N = N + 1;
    end

    b = [coefficients, fliplr(coefficients(1:end-2))];

    % Frequency response
    [H, f] = freqz(b, 1, 1024, 'whole'); % Full frequency response
    H_dB = 20*log10(abs(H));
    
    % Normalize the frequency vector
    f = f/pi;

    % Define frequency bands
    passband = f <= Fp;
    stopband = f >= Fs;

    % Calculate passband ripple and stopband attenuation
    passband_ripple = max(H_dB(passband)) - min(H_dB(passband));
    stopband_min_attenuation = min(H_dB(stopband));

    % Error between GA-generated coefficients and ideal sinc coefficients
    coeff_error = sum((b - ideal_coefficients).^2);

    % Fitness function: a combination of ripple and attenuation with penalties
    % fitness = passband_ripple + abs(stopband_min_attenuation + As) + passband_penalty + stopband_penalty ;
    fitness = 100*coeff_error;
end