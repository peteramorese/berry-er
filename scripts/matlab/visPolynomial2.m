function visPolynomial2(coeffs, ws_bounds)
    N = 100;
    [x, y] = meshgrid(linspace(ws_bounds(1), ws_bounds(2), N), linspace(ws_bounds(3), ws_bounds(4), N));

    deg = sqrt(length(coeffs));
    coeffs = reshape(coeffs, deg, deg);
    
    % Evaluate the 2D polynomial
    z = polyval2(coeffs, x, y, 1.2);
    
    % Plot the 2D polynomial
    surf(x, y, z, 'EdgeColor', 'none');
    xlabel('x_0');
    ylabel('x_1');
    zlabel('B(x)');
    % title('Visualization of 2D Polynomial');
    colormap("parula"); 
    % colorbar; % Add a colorbar to show the scale
end

