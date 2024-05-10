clear;
close all; clc;

% p_1 = [1, 2; 3, -2];
% 
% p_2 = [3, 0.5; 1, -0.5];

p_1 = [2.26806, -17.6495, 60.0409, -116.617, 141.44, -109.688, 53.1114, 13.1215, -20.6003];
p_3 = [1.56123, -2.96266, -52.0641, 340.117, -1002.13, 1761.75, -1982.68, 1416.46, -570.522, 397.115, -885.079, 885.668, -300.572]
% p_2 = [1.21212, -2.42424, 3.0303];

visPolynomial(p_1, 'red')
hold on
visPolynomial(p_3, 'green')
% visPolynomial(p_2, 'green')

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