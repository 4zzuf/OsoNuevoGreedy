function B_alpha = get_alpha_perm(n_rows, n_cols, position_array)

pairs = nchoosek(1:n_cols, 2);
n_pairs = min(n_rows, size(pairs, 1));
pairs = pairs(1:n_pairs, :);

positive = position_array(pairs(:, 1));
negative = position_array(pairs(:, 2));

B_alpha = zeros(n_rows, n_cols);
rows = (1:n_pairs).';
B_alpha(sub2ind([n_rows, n_cols], rows, positive)) = 1;
B_alpha(sub2ind([n_rows, n_cols], rows, negative)) = -1;

end
