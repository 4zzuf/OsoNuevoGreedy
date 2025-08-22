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

% Check availability of Parallel Computing Toolbox
use_parallel = license('test','Distrib_Computing_Toolbox');

for k = 1:alpha
    best_Cz = [];
    % Mask to avoid re-evaluating already selected indices
    available = true(1, n_max_comb);
    available(selected_idx) = false;

    metrics = inf(1, n_max_comb);
    Cz_cell = cell(1, n_max_comb);
    B_cell = cell(1, n_max_comb);
    sel_idx = selected_idx; % broadcast for parfor

    if use_parallel
        parfor j = 1:n_max_comb
            if ~available(j)
                continue;
            end
            candidate_idx = [sel_idx j];
            B_candidate = B_all(candidate_idx, :);
            new_B = [I_nr_r; B_candidate];
            Cz_candidate = new_B * H_r * Cx_r * H_r' * new_B' + new_B * Cn_r * new_B';
            metrics(j) = trace(Cz_candidate);
            Cz_cell{j} = Cz_candidate;
            B_cell{j} = B_candidate;
        end
    else
        for j = 1:n_max_comb
            if ~available(j)
                continue;
            end
            candidate_idx = [sel_idx j];
            B_candidate = B_all(candidate_idx, :);
            new_B = [I_nr_r; B_candidate];
            Cz_candidate = new_B * H_r * Cx_r * H_r' * new_B' + new_B * Cn_r * new_B';
            metrics(j) = trace(Cz_candidate);
            Cz_cell{j} = Cz_candidate;
            B_cell{j} = B_candidate;
        end
    end

    [~, best_idx] = min(metrics);
    selected_idx = [selected_idx best_idx];
    best_Cz = Cz_cell{best_idx};
    selected_B = B_cell{best_idx};
end

final_B = [I_nr_r; selected_B];
Cz_r = best_Cz;
K = diag(1 ./ sqrt(diag(Cz_r)));
final_W = [];
end
