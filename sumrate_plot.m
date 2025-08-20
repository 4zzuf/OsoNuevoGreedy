%SUMRATE_PLOT Compute average capacity and plot for all selection methods.
% This script assumes that the variables I_random, I_select, I_full,
% I_withoutB, capacities_exhaustive, SNR_dB and channel_realizations exist
% in the current workspace. It averages the capacities over channel
% realizations and generates a comparative plot including the exhaustive
% search results.

% Promediar las capacidades
I_random_av = sum(I_random, 2) / channel_realizations;
I_select_av = sum(I_select, 2) / channel_realizations;
I_full_av = sum(I_full, 2) / channel_realizations; % Promedio de la capacidad de la red completa
I_withoutB_av = sum(I_withoutB, 2) / channel_realizations;
I_exhaustive_av = sum(capacities_exhaustive, 2) / channel_realizations;

%%% Graficar la capacidad optimizada para ambas selecciones %%%
figure;
plot(SNR_dB, I_random_av, 'g-s', 'LineWidth', 2, 'DisplayName', 'Random Comparator Network');
hold on;
plot(SNR_dB, I_select_av, 'b-o', 'LineWidth', 2, 'DisplayName', 'Proposed Optimized Comparator Network');
plot(SNR_dB, I_full_av, 'r-x', 'LineWidth', 2, 'DisplayName', 'Full Comparator Network');
plot(SNR_dB, I_withoutB_av, 'c-x', 'LineWidth', 2, 'DisplayName', 'Without Comparator Network');
plot(SNR_dB, I_exhaustive_av, 'k-^', 'LineWidth', 2, 'DisplayName', 'Optimal (Exhaustive Search)');
xlabel('SNR (dB)', 'Interpreter', 'latex');
ylabel('Capacidad (bits/s/Hz)', 'Interpreter', 'latex');
grid on;
title('Capacidad optimizada vs. SNR', 'Interpreter', 'latex');
legend('show');

matlab2tikz('filename', 'graph201.tex');
