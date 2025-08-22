function x_qpsk = get_modulation(nt, n_bits, bits_symbol, constelation_points, gray_code_data, x)

n_symbols = n_bits / bits_symbol;

% Map Gray-coded bit patterns to constellation indices
gray_decimal = bi2de(gray_code_data, 'left-msb');
map = zeros(1, length(constelation_points));
map(gray_decimal + 1) = 1:length(constelation_points);

% Group bits into symbols and convert to indices
bit_pairs = reshape(x.', bits_symbol, []).';
idx = bi2de(bit_pairs, 'left-msb');
idx_matrix = reshape(idx, n_symbols, nt).';

x_qpsk = constelation_points(map(idx_matrix + 1));

end
