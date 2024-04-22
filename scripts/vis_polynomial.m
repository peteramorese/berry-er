clear; close all; clc;

p_1 = [1, 2; 3, -2];

p_2 = [3, 0.5; 1, -0.5];

% p_1 = [1, -1.7, 0.8, 0, 2];
% 
% p_2 = [0.6952, -0.5264, 0.2048, 0.1024, 0.0512];

visPolynomial(p_1, 'red')
hold on
visPolynomial(p_2, 'green')

function visPolynomial(coeffs, color)
    if size(coeffs, 1) == 1
        x = linspace(0, 1, 100);
        
        % Evaluate the polynomial
        y = polyval(flip(coeffs), x);
        
        % Plot the polynomial
        plot(x, y, 'b', 'LineWidth', 2, 'Color', color);
        xlabel('x');
        ylabel('y');
        title('Visualization of Polynomial');
        grid on;
    else
        [x, y] = meshgrid(linspace(0, 1, 100), linspace(0, 1, 100));
        
        % Evaluate the 2D polynomial
        z = polyval2(coeffs, x, y);
        
        % Plot the 2D polynomial
        surf(x, y, z);
        xlabel('x');
        ylabel('y');
        zlabel('z');
        title('Visualization of 2D Polynomial');
        colormap(jet); % Change the colormap if desired
        % colorbar; % Add a colorbar to show the scale
    end
end

function z = polyval2(coeffs, x, y)
    z = zeros(size(x));
    for grid_i = 1:size(z, 1)
        for grid_j = 1:size(z, 2)
            for j = 1:size(coeffs, 2)
                
                z(grid_i, grid_j) = z(grid_i, grid_j) + y(grid_i, grid_j)^(j - 1) * polyval(flip(coeffs(:, j)), x(grid_i, grid_j));
            end
        end
    end
end