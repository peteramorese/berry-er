function visPolynomial(coeffs, color)
    x = linspace(0, 1, 100);
    
    % Evaluate the polynomial
    y = polyval(flip(coeffs), x);
    
    % Plot the polynomial
    plot(x, y, 'b', 'LineWidth', 2, 'Color', color);
    xlabel('x');
    ylabel('y');
    title('Visualization of Polynomial');
    grid on;
end