%%% Parámetros %%%
Nt = 4; % 
Nr = 4; % 
SNR_dB = -10:2:20; % 
SNR = 10.^(SNR_dB / 10); % 
sigma_x = 1; % 
I_Nr = eye(Nr);
I_Nr_r = eye(2 * Nr);
I_Nt = eye(Nt);
I_Nt_r = eye(2 * Nt);
Cx = sigma_x^2 * I_Nt; 
Cx_r = (1 / 2) * sigma_x^2 * I_Nt_r; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
channel_realizations =800;
noise_realizations = 100;
alpha = 2 * Nr; % 
alpha_mse=2*Nr;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
n_max_comb = sum(1:2 * Nr - 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Precompute comparator permutations and reuse them across iterations
perm_file = 'perm_matrices.mat';
if exist(perm_file, 'file')
    load(perm_file, 'B_all_indexes', 'B_all');
else
    B_all_indexes = get_all_perm(n_max_comb, 2 * Nr);
    B_all = get_total_perm(2 * Nr);
    save(perm_file, 'B_all_indexes', 'B_all');
end

% Bits y Modulación
n_bits = 100; % Número de bits por usuario
bits_symbol = 2; % Bits por símbolo
n_symbols = n_bits / bits_symbol; % Número de símbolos
M = 2^bits_symbol; % Modulación QPSK
% Constelación QPSK
constelation_points = exp(1i * pi * (2 * (1:M) - 1) / M);
constelation_symbols = flipud((dec2bin(M - 1:-1:0) - '0'));
gray_code_data = gray_code(constelation_symbols, M);
% Transmitir bits QPSK
x = (randi(2, [Nt, n_bits]) - 1);
x_qpsk = get_modulation(Nt, n_bits, bits_symbol, constelation_points, gray_code_data, x);
x_r = [real(x_qpsk); imag(x_qpsk)];
% Inicialización del BER
BER_LRA_MMSE_rand = zeros(1, length(SNR));
BER_LRA_MMSE_full = zeros(1, length(SNR));
BER_MMSE_1bit = zeros(1, length(SNR));
BER_LMMSE = zeros(1, length(SNR));
BER_LRA_MMSE_g = zeros(1, length(SNR));
BER_LRA_MMSE_Optimized = zeros(1, length(SNR));
BER_LRA_MMSE_s = zeros(1, length(SNR));
% Bucle de simulación
for h = 1:channel_realizations
    % Canal Rayleigh
    H = (randn(Nr, Nt) + 1i * randn(Nr, Nt)) / sqrt(2);
    H_r = [real(H) -imag(H); imag(H) real(H)];
    
    %%%%%%%%%%%%% Varianza del enfoque H_r%%%%%%%%%%%%%%%%%%%%%
    new_H_r = abs(H_r).^2;
    var_H_r = sum(new_H_r, 2);
    [var_H_r_sort, position] = sort(var_H_r, 'descend');
    
    % Matrices de la red de comparadores
    B_prime = 1 / sqrt(2) * get_random_perm(alpha, 2 * Nr);
    B = [I_Nr_r; B_prime];
    
    % LRA-MMSE totalmente conectado
    full = Nr * (2 * Nr - 1); % Número de comparadores totalmente conectados
    B_alpha_f = 1 / sqrt(2) * get_alpha_perm(full, 2 * Nr, position);
    B_full = [I_Nr_r; B_alpha_f];
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Variance of H_r approach%%%%%
    new_H_r = abs(H_r).^2;
    var_H_r = sum(new_H_r,2);
    [var_H_r_sort,position] = sort(var_H_r,'descend');
    [var_H_r_sort_2,position_2] = sort(var_H_r,'ascend');
    %%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Comparator Network Matrices
    B_prime = 1/sqrt(2) * get_random_perm(alpha,2*Nr);
    B = [I_Nr_r ; B_prime];
    %%%%%%%
    % B_all_indexes and B_all are precomputed outside the loop
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    M_prime_full = 2 * Nr + full;
    B_alpha_f = 1/sqrt(2) * get_alpha_perm(full,2*Nr,position);
    B_full = [I_Nr_r ; B_alpha_f];
    %Ruido
    for i = 1:length(SNR)
        sigma_n = sqrt(sigma_x^2 ./ SNR(i));
        Cn = sigma_n^2 * eye(Nr); % Noise covariance
        Cn_r = (sigma_n^2 / 2) * eye(2 * Nr); % Noise covariance for extended space
        % Fully Connected LRA-MMSE
        Cn_r_full = B_full * Cn_r * B_full';
        Cz_full = B_full * H_r * Cx_r * H_r' * B_full' + Cn_r_full;
        K_full = diag(diag(Cz_full))^(-1/2);
        Czqx_full = sqrt(2 / pi) * K_full *B_full * H_r *Cx_r;
        Czq_full = 2 / pi * (asin(K_full * real(Cz_full) * K_full));
        W_LRA_MMSE_full = Czq_full^(-1) * Czqx_full;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Canal efectivo para red completa
        temp_CzR = B_full * (H_r * Cx_r * H_r') * B_full' + B_full * Cn_r * B_full';
        k_r_full = diag(1 ./ sqrt(diag(temp_CzR)));
        H_eff_r_q_full = sqrt(2 / pi) * k_r_full * B_full * H_r;
        
        % Lambda
        lambda = (2 / pi) * ( (pi/2 - 1) + (sigma_n^2 / 2) / ( (Nt * sigma_x^2 / 2) + (sigma_n^2 / 2) ) );
        
        %% Optimización CVX
        cvx_begin quiet
            variable Deltao(M_prime_full)
            maximize( 0.5 * log_det( eye(2 * Nt) + ((sigma_x^2 / 2) / lambda) * ...
                (H_eff_r_q_full.' * diag(Deltao) * H_eff_r_q_full) ) )
            subject to
                Deltao(1:2*Nr) == 1;
                0 <= Deltao <= 1;
                sum(Deltao) == 2*Nr + alpha;
        cvx_end
        
        % Selección de comparadores
        [~, sorted_indices] = maxk(Deltao, 2*Nr + alpha);
        vector_delta_0 = zeros(M_prime_full, 1);
        vector_delta_0(sorted_indices) = 1;
        selected_indices = find(vector_delta_0);
        B_select_alpha = B_full(selected_indices, :);
        
        %% Filtro LRA-MMSE Optimizado
        Cn_sel = B_select_alpha * Cn_r * B_select_alpha';
        Cz_sel = B_select_alpha * H_r * Cx_r * H_r' * B_select_alpha' + Cn_sel;
        K_sel = diag(1 ./ sqrt(diag(Cz_sel)));
        Czqx_sel = sqrt(2 / pi) * K_sel * B_select_alpha * H_r * Cx_r;
        Czq_sel = (2 / pi) * asin(K_sel * real(Cz_sel) * K_sel);
        W_LRA_MMSE_opt = pinv(Czq_sel) * Czqx_sel;
        
        %% Capacidad teórica para evaluación
        Cz_r = B_select_alpha * (H_r * Cx_r * H_r') * B_select_alpha' + B_select_alpha * Cn_r * B_select_alpha';
        k_r = diag(1 ./ sqrt(diag(Cz_r)));
        H_eff_r_q1 = sqrt(2 / pi) * k_r * B_select_alpha * H_r;
        C_eta_eff_r = (2 / pi) * (asin(k_r * Cz_r * k_r) - k_r * Cz_r * k_r) + ...
                      k_r * B_select_alpha * Cn_r * B_select_alpha' * k_r;
        %I_select(i_SNR, i_channel) = 0.5 * log2( det( eye(2*Nr + alpha) + ...
        %    pinv(real(C_eta_eff_r)) * ((sigma_x^2 / 2) * (H_eff_r_q1 * H_eff_r_q1') ) ) );
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Random LRA-MMSE
        Cn_rand = B * Cn_r * B';
        Cz_rand = B * H_r * Cx_r * H_r' * B' + Cn_rand;
        K_rand = diag(diag(Cz_rand))^(-1/2);
        Czqx_rand = sqrt(2 / pi) * K_rand * B * H_r*Cx_r ;
        Czq_rand = 2 / pi * (asin(K_rand * real(Cz_rand) * K_rand));
        W_LRA_MMSE_rand = Czq_rand^(-1) * Czqx_rand;
        % MMSE 1-bit Weights
        K = diag(diag(H * Cx * H' + Cn))^(-1/2);
        W_MMSE_1bit = (sqrt(pi/2) * Cx * H' * K) * ...
            (asin(transpose(K) * real(H * Cx * H' + Cn) * K) + ...
             1i * asin(transpose(K) * imag(H * Cx * H' + Cn) * K))^(-1);
        % LMMSE Weights
        W_LMMSE = (H' * H + (sigma_n^2 / sigma_x^2) * eye(Nt)) \ H';
        %MSE for Greedy     
            [W_LRA_MMSE_g, B_g, MSE_1]= greedy_search(B_all, alpha_mse, 2*Nr, I_Nr_r, Cn_r, H_r, Cx_r, n_max_comb, h, MSE_1);
        %%SINR
            [W_LRA_MMSE_s, B_s, MSE_1]= get_B_opt_2_seq_sinr(B_alpha_f, alpha, I_Nr_r, Cn_r, H_r, Cx_r,sigma_n, Nt, Nr);
       %%%%%%
        for n = 1:noise_realizations
            %% Noise and Received Signal
            n = (sigma_n * (randn(Nr, n_symbols) + 1i * randn(Nr, n_symbols))) / sqrt(2);
            y = H * x_qpsk + n;
            y_r = H_r * x_r + [real(n); imag(n)];
            z_r = B_select_alpha * y_r;
            y_1bit = sign(real(y)) + 1i * sign(imag(y));
            z_1bit = sign(real(z_r)) + 1i * sign(imag(z_r));
            %% Fully Connected LRA-MMSE Detection
            z_r_c = B_full * y_r;
            z_1bit_c = sign(real(z_r_c)) + 1i * sign(imag(z_r_c));
            x_til_r_3 = W_LRA_MMSE_full' * z_1bit_c;
            x_til_3 = x_til_r_3(1:Nt, :) + 1i * x_til_r_3(Nt+1:2*Nt, :);
            x_hat_3 = get_mapped(Nt, n_bits, bits_symbol, x_til_3);
            BER_LRA_MMSE_full(i) = BER_LRA_MMSE_full(i) + sum(sum(bitxor(x_hat_3, x)));
            %% Random LRA-MMSE Detection
            z_r_2 = B * y_r;
            z_1bit_2 = sign(real(z_r_2)) + 1i * sign(imag(z_r_2));
            x_til_r_2 = W_LRA_MMSE_rand' * z_1bit_2;
            x_til_2 = x_til_r_2(1:Nt, :) + 1i * x_til_r_2(Nt+1:2*Nt, :);
            x_hat_2 = get_mapped(Nt, n_bits, bits_symbol, x_til_2);
            BER_LRA_MMSE_rand(i) = BER_LRA_MMSE_rand(i) + sum(sum(bitxor(x_hat_2, x)));
            %% MMSE 1-bit Detection
            x_til_1bit = W_MMSE_1bit * y_1bit;
            x_hat_1bit = get_mapped(Nt, n_bits, bits_symbol, x_til_1bit);
            BER_MMSE_1bit(i) = BER_MMSE_1bit(i) + sum(sum(bitxor(x_hat_1bit, x)));
            %% LMMSE Detection
            x_til_lmmse = W_LMMSE * y;
            x_hat_lmmse = get_mapped(Nt, n_bits, bits_symbol, x_til_lmmse);
            BER_LMMSE(i) = BER_LMMSE(i) + sum(sum(bitxor(x_hat_lmmse, x)));
            
        
            %% Greedy Search
            z_r_6 =  B_g*y_r;
            z_1bit_6 = sign(real(z_r_6)) + 1i*sign(imag(z_r_6));
             x_til_r_6 = W_LRA_MMSE_g*z_1bit_6;
            x_til_LRA_6 = x_til_r_6(1:Nt,:) + 1i*x_til_r_6(Nt+1:2*Nt,:);
            
            % Mapping
            x_hat_LRA_6 = get_mapped(Nt,n_bits,bits_symbol,x_til_LRA_6);
            
            % Counting errors
            n_errors_7 = sum(sum(bitxor(x_hat_LRA_6,x)));
            BER_LRA_MMSE_g(i) = n_errors_7 + BER_LRA_MMSE_g(i);
           %% SINR
           z_r_seq = B_s*y_r;
            z_1bit_seq = sign(real(z_r_seq)) + 1i*sign(imag(z_r_seq));
            x_til_r_seq = W_LRA_MMSE_s.'*z_1bit_seq;
            x_til_LRA_seq = x_til_r_seq(1:Nt,:) + 1i*x_til_r_seq(Nt+1:2*Nt,:);
            
            % Mapping
            x_hat_LRA_seq = get_mapped(Nt,n_bits,bits_symbol,x_til_LRA_seq);
            
            % Counting errors
            n_errors_s = sum(sum(bitxor(x_hat_LRA_seq,x)));
            BER_LRA_MMSE_s(i) = n_errors_s + BER_LRA_MMSE_s(i);
             % Detección LRA-MMSE optimizada
        x_tilde_r = W_LRA_MMSE_opt' * z_1bit;
        x_tilde = x_tilde_r(1:Nt,:) + 1i * x_tilde_r(Nt+1:end,:);
        %%Mapping
        x_hat_LRA = get_mapped(Nt, n_bits, bits_symbol, x_tilde);
        
        % Counting Errors
        n_errors_s = sum(sum(bitxor(x_hat_LRA, x)));
        BER_LRA_MMSE_Optimized(i) = BER_LRA_MMSE_Optimized(i) + n_errors_s;
            
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:length(SNR)
        sigma_n = sqrt(sigma_x^2 / SNR(i));
        Cn_r = (sigma_n^2 / 2) * eye(2 * Nr); % Covarianza de ruido en espacio extendido
        
        disp(['Calculando para SNR = ', num2str(SNR_dB(i)), ' dB']);
        %
        Cn_r_full = B_full * Cn_r * B_full';
        Cz_full = B_full * H_r * Cx_r * H_r' * B_full' + Cn_r_full;
        K_full = diag(diag(Cz_full))^(-1/2);
        Czqx_full = sqrt(2 / pi) * K_full * B_full * H_r * Cx_r;
        Czq_full = 2 / pi * (asin(K_full * real(Cz_full) * K_full));
        W_LRA_MMSE_full = Czq_full^(-1) * Czqx_full;
        %%%%%%%%%%%%%%%%%%%%%%
        temp_CzR = B_full * (H_r * Cx_r * H_r') * B_full' + B_full * Cn_r * B_full';
        k_r_full = diag(1 ./ sqrt(diag(temp_CzR)));
        H_eff_r_q_full = sqrt(2 / pi) * k_r_full * B_full * H_r;
        
        % 
        k = alpha;
        combinations = nchoosek(1:full, k);
        max_capacity = -inf;
        min_BER = Inf;
        
        % Búsqueda exhaustiva para los mejores comparadores
        for j = 1:size(combinations, 1)
            selected_indices_exh = combinations(j, :);
            B_select_exh_prime = B_alpha_f(selected_indices_exh, :);
            B_select_exh = [I_Nr_r; B_select_exh_prime];
            
            % Capacidad y Detección para la combinación seleccionada
            Cz_r_exh = B_select_exh * (H_r * Cx_r * H_r') * B_select_exh' + B_select_exh * Cn_r * B_select_exh';
            k_r_exh = diag(1 ./ sqrt(diag(Cz_r_exh))); 
            H_eff_r_q_exh = sqrt(2 / pi) * k_r_exh * B_select_exh * H_r;           
            C_eta_eff_r_exh = (2 / pi) * (asin(k_r_exh * Cz_r_exh * k_r_exh) - k_r_exh * Cz_r_exh * k_r_exh) + ...
                               k_r_exh * B_select_exh * Cn_r * B_select_exh' * k_r_exh;
            capacity = 1 / 2 * log2(det(eye(2 * Nr + alpha) + pinv(real(C_eta_eff_r_exh)) * ...
                ((sigma_x^2 / 2) * (H_eff_r_q_exh * H_eff_r_q_exh'))));
            if capacity > max_capacity
                max_capacity = capacity;
            end
            
            % Detección y BER
            errors_total = 0;
            for n_sim = 1:10
                x = (randi(2, [Nt, n_bits]) - 1);
                x_qpsk = get_modulation(Nt, n_bits, bits_symbol, constelation_points, gray_code_data, x);
                x_r = [real(x_qpsk); imag(x_qpsk)];
                n_noise = (sigma_n * (randn(Nr, n_symbols) + 1i * randn(Nr, n_symbols))) / sqrt(2);
                y_r = H_r * x_r + [real(n_noise); imag(n_noise)];
                % Detección LRA-MMSE
                Cn_exh = B_select_exh * Cn_r * B_select_exh';
                Cz_exh = B_select_exh * H_r * Cx_r * H_r' * B_select_exh' + Cn_exh;
                K_exh = diag(1 ./ sqrt(diag(Cz_exh)));
                Czqx_exh = sqrt(2 / pi) * K_exh * B_select_exh * H_r * Cx_r;
                Czq_exh = (2 / pi) * asin(K_exh * real(Cz_exh) * K_exh);
                W_LRA_MMSE_exh = pinv(Czq_exh) * Czqx_exh;
                % Detección
                z_r_exh = B_select_exh * y_r;
                z_1bit_exh = sign(real(z_r_exh)) + 1i * sign(imag(z_r_exh));
                x_til_r_exh = W_LRA_MMSE_exh' * z_1bit_exh;
                x_til_exh = x_til_r_exh(1:Nt, :) + 1i * x_til_r_exh(Nt + 1:2 * Nt, :);
                x_hat_exh = get_mapped(Nt, n_bits, bits_symbol, x_til_exh);
                errors_total = errors_total + sum(sum(bitxor(x_hat_exh, x)));
            end
            BER_exh_candidate = errors_total / (10 * Nt * n_bits);
            if BER_exh_candidate < min_BER
                min_BER = BER_exh_candidate;
            end
        end
        
        capacities_exhaustive(i, h) = max_capacity;
        BER_exhaustive(i, h) = min_BER;
    end
end
% Promedio sobre realizaciones del canal
BER_LRA_MMSE_rand = BER_LRA_MMSE_rand / (noise_realizations * Nt * channel_realizations * n_bits);
BER_LRA_MMSE_full = BER_LRA_MMSE_full / (noise_realizations * Nt * channel_realizations * n_bits);
BER_MMSE_1bit = BER_MMSE_1bit / (noise_realizations * Nt * channel_realizations * n_bits);
BER_LMMSE = BER_LMMSE / (noise_realizations * Nt * channel_realizations * n_bits);
BER_LRA_MMSE_g = BER_LRA_MMSE_g / (noise_realizations * Nt * channel_realizations * n_bits);
BER_LRA_MMSE_Optimized = BER_LRA_MMSE_Optimized / (noise_realizations * Nt * channel_realizations * n_bits);
BER_LRA_MMSE_s = BER_LRA_MMSE_s / (noise_realizations * Nt * channel_realizations * n_bits);
I_exhaustive_av = sum(capacities_exhaustive, 2) / channel_realizations;
BER_exhaustive_av = sum(BER_exhaustive, 2) / channel_realizations;
% Graficar BER vs. SNR
figure;
semilogy(SNR_dB, BER_LRA_MMSE_rand, 'r-o', 'LineWidth', 1.5); % 
hold on;
semilogy(SNR_dB, BER_LRA_MMSE_full, 'b-', 'LineWidth', 1.5); % 
semilogy(SNR_dB, BER_MMSE_1bit, 'm-.', 'LineWidth', 1.5); % 
semilogy(SNR_dB, BER_LMMSE, 'g--', 'LineWidth', 1.5); % 
semilogy(SNR_dB, BER_LRA_MMSE_g, 'k:', 'LineWidth', 1.5); % 
semilogy(SNR_dB, BER_LRA_MMSE_Optimized, 'b-o', 'LineWidth', 1.5); % 
semilogy(SNR_dB, BER_LRA_MMSE_s, 'c-o', 'LineWidth', 1.5); % 
semilogy(SNR_dB, BER_exhaustive_av, 'm-^', 'LineWidth', 2); % 
xlabel('SNR (dB)', 'Interpreter', 'latex');
ylabel('BER', 'Interpreter', 'latex');
grid on;
title('BER vs. SNR para diferentes detectores');
legend({'LRA-MMSE (Rand)', 'LRA-MMSE (Full)', 'MMSE 1bit', 'LMMSE', 'LRA-MMSE (g)', 'LRA-MMSE (Optimized)', 'LRA-MMSE (s)', 'Búsqueda Exhaustiva LRA-MMSE'}, 'Location', 'best');
hold off;
matlab2tikz('filename', 'graphdet2.tex');