function z = polyval2(coeffs, x, y, cap)
    z = zeros(size(x));
    for grid_i = 1:size(z, 1)
        for grid_j = 1:size(z, 2)
            for j = 1:size(coeffs, 2)
                z(grid_i, grid_j) = z(grid_i, grid_j) + y(grid_i, grid_j)^(j - 1) * polyval(flip(coeffs(:, j)), x(grid_i, grid_j));
                % z(grid_i, grid_j) = min(z(grid_i, grid_j), cap);
            end
        end
    end
    z = min(z, cap);
    fprintf("Min grid value: %.5f\n", min(z,[],"all"));
end

