clc; clear; close all;
%rng(2);

nt = 4;
nr = 4;
SNR_dB = -20:5:30;
SNR = 10.^(SNR_dB/10);
Etx = 1;
I_nr = eye(nr);
I_nr_r = eye(2*nr);
I_nt = eye(nt);
Cx = I_nt;
I_nt_r = eye(2*nt);
Cx_r = (1/2)*I_nt_r;
I_nt_nr = eye(nt*nr);
I_2nt_nr = eye(2*nt*nr);
channel_realizations = 10;
noise_realizations = 10;
alpha = 2*nr;
alpha_mse = 2*nr;
n_max_comb = 0;
for a = 1:2*nr-1
    n_max_comb = n_max_comb + a;
end
alpha_fully = n_max_comb;
n_bits = 100;
bits_symbol = 2;
n_symbols = n_bits/bits_symbol;
M = 2^bits_symbol;
constelation_points = zeros(1,M);
for t = 1:M
    constelation_points(t) = exp(1i*pi*(2*t-1)/M);
end
constelation_symbols = flipud((dec2bin(M-1:-1:0)-'0'));
gray_code_data = gray_code(constelation_symbols, M);

pilot_length = 2*nt;
o = orth(randn(pilot_length,pilot_length) + 1i*randn(pilot_length,pilot_length));
p = o(1:nt,:)*sqrt(pilot_length);
p_tilde = kron(p.',I_nr);
p_tilde_r = [real(p_tilde) -imag(p_tilde); imag(p_tilde) real(p_tilde)];
new = p_tilde_r*p_tilde_r.';
I_tau = eye(pilot_length);
I_tau_nr = eye(pilot_length*nr);
I_tau_2nr = eye(2*pilot_length*nr);
x = (randi(2,[nt, n_bits])-1);
x_qpsk = get_modulation(nt, n_bits, bits_symbol, constelation_points, gray_code_data, x);
x_r = [real(x_qpsk); imag(x_qpsk)];

BER_LRA_MMSE = zeros(length(nr), length(SNR));
BER_LMMSE = zeros(length(nr), length(SNR));
BER_LRA_MMSE_2 = zeros(length(nr), length(SNR));
BER_LMMSE_2 = zeros(length(nr), length(SNR));
BER_LRA_MMSE_3 = zeros(length(nr), length(SNR));
BER_LMMSE_3 = zeros(length(nr), length(SNR));
BER_LRA_MMSE_4 = zeros(length(nr), length(SNR));
BER_LMMSE_4 = zeros(length(nr), length(SNR));

MSE_data = zeros(1, length(SNR));
MSE_data_2 = zeros(1, length(SNR));
MSE_data_3 = zeros(1, length(SNR));
MSE_data_4 = zeros(1, length(SNR));

n_MSE_data = zeros(channel_realizations, noise_realizations, length(SNR));
n_MSE_data_2 = zeros(channel_realizations, noise_realizations, length(SNR));
n_MSE_data_3 = zeros(channel_realizations, noise_realizations, length(SNR));
n_MSE_data_4 = zeros(channel_realizations, noise_realizations, length(SNR));

n_MSE_mean = zeros(1, length(SNR));
n_MSE_mean_2 = zeros(1, length(SNR));
n_MSE_mean_3 = zeros(1, length(SNR));
n_MSE_mean_4 = zeros(1, length(SNR));

r2_tmp = zeros(2*nt, 1);
r2_tmp_2 = zeros(2*nt, 1);
r2_tmp_3 = zeros(2*nt, 1);
r2_tmp_4 = zeros(2*nt, 1);
r3_tmp = zeros(2*nt, 1);

Ir_k = zeros(2*nt, noise_realizations, channel_realizations, length(SNR));
Ir_k_2 = zeros(2*nt, noise_realizations, channel_realizations, length(SNR));
Ir_k_3 = zeros(2*nt, noise_realizations, channel_realizations, length(SNR));
Ir_k_4 = zeros(2*nt, noise_realizations, channel_realizations, length(SNR));

Ir_mean = zeros(2*nt, length(SNR));
Ir_final = zeros(1, length(SNR));

n_dif_channels = 100;
for m = 1:n_dif_channels
    H(:,:,m) = (randn(nr, nt) + 1i*randn(nr, nt))/sqrt(2);
    H_r(:,:,m) = [real(H(:,:,m)) -imag(H(:,:,m)); imag(H(:,:,m)) real(H(:,:,m))];
end

H_r_this = H_r;

for h = 1:channel_realizations
    vec_h = reshape(H(1:nr,:,h), nt*nr, 1);
    Rh = I_nt_nr;
    vec_h_r = [real(vec_h); imag(vec_h)];
    Rh_2 = (1/2)*I_2nt_nr;

    % B matrix (Greedy Search)
    B_all_indexes = get_all_perm(n_max_comb, 2*nr);
    B_all = get_total_perm(2*nr);

    % Beff matrix (Greedy Search)
    B_r_prime_3 = B_all(:, 1:nr);
    Beff_r_prime_3 = kron(B_r_prime_3, I_tau);
    B_i_prime_3 = B_all(:, nr+1:2*nr);
    Beff_i_prime_3 = kron(B_i_prime_3, I_tau);
    Beff_prime_3 = [Beff_r_prime_3 Beff_i_prime_3];
    Beff_3 = [I_tau_2nr; Beff_prime_3];

    for i = 1:length(SNR)
        sigma_w = sqrt(Etx./SNR(i));
        Cw = sigma_w^2*eye(nr);
        Cw_new = ((sigma_w^2)/2)*eye(2*nr);

        Cw_r_3 = B_all*Cw_new*B_all.';
        Cn_p = (sigma_w^2)*I_tau_nr;
        Cn_r_p = ((sigma_w^2)/2)*I_tau_2nr;

        % Achievable Rate (Greedy)
        A_d_4 = sqrt(2/pi)*B_all;
        Cnq_d_4 = Cw_new - A_d_4*Cw_r_3*A_d_4;
        for n = 1:noise_realizations
            w = (sigma_w*(randn(nr, pilot_length + n_symbols) + 1i*randn(nr, pilot_length + n_symbols)))/sqrt(2);
            w_r = [real(w(:,pilot_length+1:pilot_length+n_symbols)); imag(w(:,pilot_length+1:pilot_length+n_symbols))];
            y_r = H_r_this(:,:,h)*x_r + w_r;
            z_1bit_4 = sign(y_r);

            z_r_4 = B_all*y_r;
            z_1bit_5 = sign(z_r_4);
            
            n_MSE_data_4(n, h, i) = (norm((z_1bit_5 - vec_h), 2))^2;
        end
    end
end

% Save results
save('resultado.mat');
disp(h);
end

BER_LRA_MMSE = BER_LRA_MMSE/(noise_realizations*nt*channel_realizations*n_bits);
BER_LMMSE = BER_LMMSE/(noise_realizations*nt*channel_realizations*n_bits);
BER_LRA_MMSE_2 = BER_LRA_MMSE_2/(noise_realizations*nt*channel_realizations*n_bits);
BER_LMMSE_2 = BER_LMMSE_2/(noise_realizations*nt*channel_realizations*n_bits);
% Plotting BER vs SNR
figure(1)
semilogy(SNR_dB, BER_LRA_MMSE(1,:),'g-');
hold on
semilogy(SNR_dB, BER_LMMSE(1,:),'g--');
hold off
xlabel('SNR (dB)');
ylabel('BER');
grid on;
title('BER vs. SNR');
legend('LRA-MMSE', 'LMMSE');
