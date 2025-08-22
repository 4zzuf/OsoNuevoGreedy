
function m_comb = get_all_perm(n_max_comb, cols)
%GET_ALL_PERM Generate up to the first N unique pair combinations
%   m_comb = GET_ALL_PERM(cols) returns all unique combinations of two
%   indices from 1 to cols.
%   m_comb = GET_ALL_PERM(n_max_comb, cols) returns only the first
%   n_max_comb unique combinations of two indices from 1 to cols.

% Allow calling with a single argument that specifies only the number of
% columns. In that case, generate all possible combinations.
if nargin == 1
    cols = n_max_comb;
    n_max_comb = cols * (cols - 1) / 2;
end

% Total number of unique pairs and validation of the requested amount
total = cols * (cols - 1) / 2;
if n_max_comb > total
    error('n_max_comb exceeds the total number of combinations (%d).', total);
end

% Preallocate and fill the combinations using nested loops
m_comb = zeros(n_max_comb, 2);
row = 1;
for a = 1:cols-1
    for b = a + 1:cols
        m_comb(row, :) = [a, b];
        row = row + 1;
        if row > n_max_comb
            return; % Stop once the desired number of combinations is reached
        end
    end
end

end
