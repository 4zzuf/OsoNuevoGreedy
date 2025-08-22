function B_alpha = get_alpha_perm(n_rows, n_cols, position_array) 

B_alpha = zeros([n_rows, n_cols]);
        
B_alpha_row = 1;
for step=1:n_cols-1
    cont = 1;
    while cont < n_cols-step+1
        positive_part = position_array(cont);
        negative_part = position_array(cont+step);
        B_alpha(B_alpha_row,positive_part) = 1;
        B_alpha(B_alpha_row,negative_part) = -1;

        if B_alpha_row == n_rows
            return;
        end

        B_alpha_row = B_alpha_row+1;
        cont = cont + 1;

    end
end
end