function [B_es, K_es, Cz_es] = es_search(B_all, alpha, I_Nr_r, Cn_r, H_r, Cx_r)
%ES_SEARCH Exhaustive comparator network search.
%   [B_es, K_es, Cz_es] = ES_SEARCH(B_all, alpha, I_Nr_r, Cn_r, H_r, Cx_r)
%   performs an exhaustive search over all combinations of ALPHA rows from
%   the candidate comparator matrix B_all. The combination that maximizes
%   the mutual information is returned as B_es together with its
%   normalisation matrix K_es and the corresponding receive covariance
%   Cz_es.
%
%   This approach has factorial complexity and is intended only for small
%   problem sizes.

% number of available comparator rows
n_candidates = size(B_all, 1);

% all combinations of alpha rows
comb_indices = nchoosek(1:n_candidates, alpha);

% recover sigma_x^2 from Cx_r = (1/2)*sigma_x^2*I
sigma_x_sq = 2 * Cx_r(1,1);

best_capacity = -inf;
B_es = [];
K_es = [];
Cz_es = [];

for idx = 1:size(comb_indices, 1)
    sel_rows = comb_indices(idx, :);
    B_candidate = [I_Nr_r; B_all(sel_rows, :)];

    Cz_candidate = B_candidate * H_r * Cx_r * H_r' * B_candidate' + ...
                   B_candidate * Cn_r * B_candidate';
    K_candidate = diag(1 ./ sqrt(diag(Cz_candidate)));

    H_eff = sqrt(2/pi) * K_candidate * B_candidate * H_r;
    C_eta = (2/pi) * (asin(K_candidate * Cz_candidate * K_candidate) - ...
                     K_candidate * Cz_candidate * K_candidate) + ...
            K_candidate * B_candidate * Cn_r * B_candidate' * K_candidate;

    m = size(B_candidate, 1);
    capacity = 0.5 * log2(det(eye(m) + pinv(real(C_eta)) * ...
                        ((sigma_x_sq/2) * H_eff * H_eff')));

    if capacity > best_capacity
        best_capacity = capacity;
        B_es = B_candidate;
        K_es = K_candidate;
        Cz_es = Cz_candidate;
    end
end

end
