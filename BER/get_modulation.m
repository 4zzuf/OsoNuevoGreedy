function x_qpsk = get_modulation(nt,n_bits,bits_symbol,constelation_points,gray_code_data,x)

x_qpsk = zeros([nt, n_bits/bits_symbol]);

for row=1:nt
    col = 1;
    for counter=1:bits_symbol:n_bits
        for counter2=1:length(constelation_points)
            
            b1 = x(row,counter);
            b2 = x(row,counter+1);
            symbol = [b1,b2];
            
            if hamming_distance(symbol,gray_code_data(counter2,:))==0
                x_qpsk(row,col) = constelation_points(counter2);
                col = col + 1;
            end
            
        end
    end
end

end