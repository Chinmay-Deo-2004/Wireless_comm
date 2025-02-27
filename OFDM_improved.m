close all; clc;

%% OFDM Simulation
% Initializing parameters
Nsc = input('OFDM symbol size (Number of subcarriers) N = '); 
M = input('Modulation order M = '); 
Nsmb = input('Number of OFDM symbols to be simulated = '); 
Ne = 3000; % Maximum error bits
step_size = input('SNR step size [dB] = '); 
SNR_start = 0; 
SNR_end = input('Last value of SNR [dB] = ');

% Derived parameters
bits_per_symbol = log2(M); % Number of bits per QAM symbol
disp(['No. of bits per symbol = ', num2str(bits_per_symbol)]);

% SNR range
snr_range = SNR_start:step_size:SNR_end;

% Convolutional code parameters
trellis = poly2trellis(7, [133 171]); % Rate 1/2 convolutional code
code_rate = 1 / 2;

% Placeholder for BER results
ber_results = zeros(1, length(snr_range));

%% Monte Carlo simulation
for idx = 1:length(snr_range)
    snr = snr_range(idx);
    disp(['Processing SNR = ', num2str(snr), ' dB']);

    % Initialize error and symbol counters
    num_errors = 0;
    num_symbols = 0;

    while num_errors < Ne && num_symbols < Nsmb
        % Transmitter
        data_bits = randi([0, M-1], 1, Nsc); % Random data generation
        coded_bits = convenc(data_bits, trellis); % Convolutional encoding
        
        modulated_signal = 2 * data - 1; % BPSK modulation (Map 0 -> -1, 1 -> +1)

        %modulated_signal = qammod(data, M, 'UnitAveragePower', true); % QAM modulation

        ofdm_signal = ifft(modulated_signal, Nsc); % OFDM modulation (IFFT)

        % Rayleigh fading channel
        h = (randn(1, Nsc) + 1j * randn(1, Nsc)) / sqrt(2); % Rayleigh fading
        faded_signal = h .* ofdm_signal; % Apply fading

        % Add noise (AWGN)
        noisy_signal = awgn(faded_signal, snr, 'measured');

        % Receiver
        % Equalize the received signal (remove channel effect)
        h_est = mean(h, 2); % Channel estimation (simplified)
        received_signal = noisy_signal ./ h; 

        % OFDM demodulation (FFT)
        demodulated_signal = fft(received_signal, Nsc);

        % Demodulate symbols
        demodulated_data = real(demodulated_signal) > 0; % BPSK demodulation

        % Viterbi decoding
        decoded_bits = vitdec(received_bits, trellis, 34, 'trunc', 'hard'); % Decode with Viterbi

        % Bit error calculation
        [bit_errors, ~] = biterr(data, demodulated_data);
        num_errors = num_errors + bit_errors;
        num_symbols = num_symbols + 1;
    end

    % Calculate BER for current SNR
    ber_results(idx) = num_errors / (bits_per_symbol * num_symbols);
end

%% Plot Results
figure;
semilogy(snr_range, ber_results, '-ok', 'LineWidth', 2, ...
         'MarkerFaceColor', 'k', 'MarkerSize', 8, 'MarkerEdgeColor', 'k');
grid on;
title('OFDM Bit Error Rate vs SNR (Rayleigh Fading)');
xlabel('SNR [dB]');
ylabel('Bit Error Rate (BER)');
legend(['BER, N = ', num2str(Nsc), ', ', num2str(M), '-QAM']);