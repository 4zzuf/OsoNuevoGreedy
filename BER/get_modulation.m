function x_qpsk = get_modulation(nt, n_bits, bits_symbol, constelation_points, gray_code_data, x)

n_symbols = n_bits / bits_symbol;

% Map Gray-coded bit patterns to constellation indices
gray_decimal = binary_to_decimal(gray_code_data);
map = zeros(1, length(constelation_points));
map(gray_decimal + 1) = 1:length(constelation_points);

% Group bits into symbols and convert to indices
bit_pairs = reshape(x.', bits_symbol, []).';
idx = binary_to_decimal(bit_pairs);
idx_matrix = reshape(idx, n_symbols, nt).';

x_qpsk = constelation_points(map(idx_matrix + 1));

end

function dec = binary_to_decimal(bits)
% Convert rows of a binary matrix to decimal values assuming the left-most
% column is the most significant bit. This function replaces the
% Communications Toolbox `bi2de` dependency so that the code runs in
% environments where that toolbox is unavailable.

    weights = 2 .^ (size(bits, 2) - 1:-1:0);
    dec = bits * weights.';
end
