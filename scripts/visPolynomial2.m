function visPolynomial2(coeffs, ws_bounds)
    [x, y] = meshgrid(linspace(ws_bounds(1), ws_bounds(2), 100), linspace(ws_bounds(3), ws_bounds(4), 100));

    deg = sqrt(length(coeffs));
    coeffs = reshape(coeffs, deg, deg);
    
    % Evaluate the 2D polynomial
    z = polyval2(coeffs, x, y, 10);
    
    % Plot the 2D polynomial
    surf(x, y, z);
    xlabel('x');
    ylabel('y');
    zlabel('z');
    title('Visualization of 2D Polynomial');
    colormap(jet); 
    % colorbar; % Add a colorbar to show the scale
end

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
end