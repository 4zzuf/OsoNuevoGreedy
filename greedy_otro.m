function [final_W, final_B, K, Cz_r] = greedy_otro(B_all, alpha, I_nr_r, Cn_r, H_r, Cx_r, n_max_comb)
%GREEDY_OTRO Selects comparator network rows using a simple greedy search.
%   This function implements a basic greedy algorithm that incrementally
%   builds the comparator matrix by selecting the rows from B_all that
%   minimise the trace of the resulting received covariance matrix.
%
%   [final_W, final_B, K, Cz_r] = GREEDY_OTRO(B_all, alpha, I_nr_r, Cn_r,
%   H_r, Cx_r, n_max_comb) returns the selected weight matrix final_B, the
%   normalisation matrix K and the receive covariance Cz_r. final_W is left
%   empty since it is not required by untitled1.m but kept for interface
%   compatibility with GREEDY_SEARCH.
%
%   Inputs:
%       B_all      - Matrix with all candidate comparator rows.
%       alpha      - Number of rows to be selected from B_all.
%       I_nr_r     - Identity matrix of size 2*Nr.
%       Cn_r       - Noise covariance matrix in the real domain.
%       H_r        - Real-valued channel matrix.
%       Cx_r       - Transmit covariance matrix in the real domain.
%       n_max_comb - Total number of candidate rows available in B_all.
%
%   Outputs:
%       final_W    - Empty matrix (kept for interface compatibility).
%       final_B    - Greedily selected comparator matrix including the
%                    identity block [I_nr_r; B_selected].
%       K          - Diagonal normalisation matrix based on Cz_r.
%       Cz_r       - Receive covariance matrix for the selected structure.
%
%   The approach is intentionally lightweight and focuses on demonstrating
%   how the new greedy routine can be integrated into the main simulation.

% Initialise set of selected indices
selected_idx = [];
selected_B = [];

for k = 1:alpha
    best_idx = 0;
    best_metric = inf;
    best_Cz = [];
    % Evaluate each unselected row
    for j = 1:n_max_comb
        if ismember(j, selected_idx)
            continue;
        end
        candidate_idx = [selected_idx j];
        B_candidate = B_all(candidate_idx, :);
        new_B = [I_nr_r; B_candidate];
        Cz_candidate = new_B * H_r * Cx_r * H_r' * new_B' + new_B * Cn_r * new_B';
        metric = trace(Cz_candidate);
        if metric < best_metric
            best_metric = metric;
            best_idx = j;
            best_Cz = Cz_candidate;
            selected_B = B_candidate;
        end
    end
    selected_idx = [selected_idx best_idx];
end

final_B = [I_nr_r; selected_B];
Cz_r = best_Cz;
K = diag(1 ./ sqrt(diag(Cz_r)));
final_W = [];
end
