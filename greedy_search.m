function [final_W, final_B, K, Cz_r] = greedy_search(B_all, alpha, I_nr_r, Cn_r, H_r, Cx_r, n_max_comb)
% Alpha first rows of B_all
B_alpha = B_all(1:alpha,:);
% Precompute constant matrices
HCxH = H_r*Cx_r*H_r';

% Computing the MSE for the first alpha rows
new_B = [I_nr_r ; B_alpha];
Cn_rand = new_B*Cn_r*new_B';
K = diag(diag(new_B*HCxH*new_B'+Cn_rand).^(-1/2));
%Czqx = sqrt(pi/2)*Cx_r*H_r'*new_B'*K;
Czqx = sqrt(pi/2)*(1/2)*Cx_r*H_r'*new_B'*K;

% Compute Czq ensuring numerical stability
Czq_matrix = asin(transpose(K) * real(new_B * HCxH * new_B' + Cn_rand) * K) ...
    + 1i * asin(transpose(K) * imag(new_B * HCxH * new_B' + Cn_rand) * K);

if size(Czq_matrix, 1) == size(Czq_matrix, 2) && cond(Czq_matrix) < 1e10
    Czq = inv(Czq_matrix);
else
    Czq = pinv(Czq_matrix);
end

W = Czqx * Czq;
%equivalent objective MSE
MSE_data = -2*real(trace(Czqx*W')) + trace(W*Czq*W');
%total MSE
MSE_data_2 = trace(Cx_r) - 2*real(trace(Czqx*W')) + trace(W*Czq*W');
% Saving the MSE with the first alpha rows
lower_value = MSE_data;
final_W = W;
final_B = new_B;
lower_value_2 = MSE_data_2;

% Saving the indexes that are already in B_alpha
opt_index_rows = zeros(1,n_max_comb);
opt_index_rows(1,1:alpha) = 1;
% Track the indices from B_all that compose B_alpha
selected_indices = 1:alpha;

% Check availability of Parallel Computing Toolbox
use_parallel = license('test','Distrib_Computing_Toolbox');

for i=1:alpha
    % Getting the rows that are not going to be permuted and freeze them
    B_frozen = B_alpha;
    B_frozen(i,:) = [];

    % Preallocate containers for parfor/for loop
    MSE_values  = inf(1, n_max_comb);
    MSE2_values = inf(1, n_max_comb);
    W_cell      = cell(1, n_max_comb);
    B_cell      = cell(1, n_max_comb);
    Balpha_cell = cell(1, n_max_comb);

    % Local copy for broadcast inside parfor
    opt_mask = opt_index_rows;

    if use_parallel
        parfor j=1:n_max_comb
            if opt_mask(j) == 1
                continue;
            end
            B_perm_row = B_all(j,:);
            B_result = [B_frozen(1:i-1,:); B_perm_row; B_frozen(i:end,:)];

            % Computing the MSE value
            new_B = [I_nr_r ; B_result];
            Cw_r = new_B*Cn_r*new_B';

            K = diag(diag(new_B*HCxH*new_B'+Cw_r).^(-1/2));
            Czqx = sqrt(pi/2)*(1/2)*Cx_r*H_r'*new_B'*K;
            Czq = pinv(asin(transpose(K)*real(new_B*HCxH*new_B'+Cw_r)*K) + 1i*asin(transpose(K)*imag(new_B*HCxH*new_B'+Cw_r)*K));
            W = Czqx*Czq;

            MSE_values(j)  = -2*real(trace(Czqx*W')) + trace(W*Czq*W');
            MSE2_values(j) = trace(Cx_r) - 2*real(trace(Czqx*W')) + trace(W*Czq*W');
            W_cell{j}      = W;
            B_cell{j}      = new_B;
            Balpha_cell{j} = B_result;
        end
    else
        for j=1:n_max_comb
            if opt_mask(j) == 1
                continue;
            end
            B_perm_row = B_all(j,:);
            B_result = [B_frozen(1:i-1,:); B_perm_row; B_frozen(i:end,:)];

            % Computing the MSE value
            new_B = [I_nr_r ; B_result];
            Cw_r = new_B*Cn_r*new_B';

            K = diag(diag(new_B*HCxH*new_B'+Cw_r).^(-1/2));
            Czqx = sqrt(pi/2)*(1/2)*Cx_r*H_r'*new_B'*K;
            Czq = pinv(asin(transpose(K)*real(new_B*HCxH*new_B'+Cw_r)*K) + 1i*asin(transpose(K)*imag(new_B*HCxH*new_B'+Cw_r)*K));
            W = Czqx*Czq;

            MSE_values(j)  = -2*real(trace(Czqx*W')) + trace(W*Czq*W');
            MSE2_values(j) = trace(Cx_r) - 2*real(trace(Czqx*W')) + trace(W*Czq*W');
            W_cell{j}      = W;
            B_cell{j}      = new_B;
            Balpha_cell{j} = B_result;
        end
    end

    % Determine best candidate for current row
    [current_min, best_j] = min(MSE_values);
    if current_min < lower_value
        lower_value = current_min;
        final_W = W_cell{best_j};
        final_B = B_cell{best_j};
        B_alpha = Balpha_cell{best_j};
        opt_index_rows(selected_indices(i)) = 0;
        selected_indices(i) = best_j;
        opt_index_rows(best_j) = 1;
    end
    lower_value_2 = min(lower_value_2, min(MSE2_values));
end

% Final Cz_r and K
Cz_r = final_B*HCxH*final_B' + final_B*Cn_r*final_B';
K = diag(diag(Cz_r).^(-1/2));

end
