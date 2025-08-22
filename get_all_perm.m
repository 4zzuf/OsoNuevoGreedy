function m_comb = get_all_perm(n_max_comb, cols)
%GET_ALL_PERM Generate all unique pair combinations
%   m_comb = GET_ALL_PERM(cols) returns all unique combinations
%   of two indices from 1 to cols.
%   m_comb = GET_ALL_PERM(n_max_comb, cols) returns up to n_max_comb
%   unique combinations of two indices from 1 to cols.

% Handle single-argument call by treating it as cols
if nargin == 1
    cols = n_max_comb;
    n_max_comb = [];
end

% Generate all combinations of two indices without explicit loops
m_comb = nchoosek(1:cols, 2);

% Truncate if a maximum number of combinations is specified
if nargin == 2
    total_combs = size(m_comb, 1);
    if n_max_comb < total_combs
        m_comb = m_comb(1:n_max_comb, :);
    elseif n_max_comb > total_combs
        warning('n_max_comb exceeds the total number of combinations (%d). Returning all combinations.', total_combs);
    end
end
end
