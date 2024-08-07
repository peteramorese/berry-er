clear; close all;

file_directory = "../../build/bin/";

N = 10;
d = 2;
% ws_bounds = [0, 1, 0, 1];
ws_bounds = [-1.2, .7, -.7, .7];



A = readMatrixFromFile(file_directory + "A.txt");
b = readMatrixFromFile(file_directory + "b.txt");
% Phi_inv = readMatrixFromFile(file_directory + "Phi_inv.txt");

b = [b; 0; 0];
gamma_cond = zeros(1, size(A,2));
gamma_cond(end) = 1;
eta_cond = zeros(1, size(A,2));
eta_cond(end-1) = 1;

A = [A; gamma_cond; eta_cond];

objective_vec = zeros(size(A,2), 1);
objective_vec(end - 1) = 1.0;
objective_vec(end) = N;
tic

% % options = optimoptions('linprog', 'Algorithm', 'interior-point');
options = optimoptions('linprog', ...
    'Algorithm', 'interior-point', ...
    'OptimalityTolerance', 1e-4, ... % Optimality tolerance
    'ConstraintTolerance', 1e-8);%, ... % Constraint tolerance (feasibility tolerance)
    % 'Display', 'iter');              % Display iteration information

vars = linprog(objective_vec, -A, -b, [], [], [], [], options);

% vars = linprog(objective_vec, -A, -b);

t = toc;
fprintf("Time: %.3f\n", t)
coeffs = vars(1:end-2);
eta = vars(end - 1);
gamma = vars(end);
fprintf("Eta: %.3f\n", eta)
fprintf("Gamma: %.3f\n", gamma)
fprintf("Probability of safety: %.3f\n", 1 - (eta + N * gamma))

if d == 1
    figure;
    visPolynomial(coeffs, "red")
elseif d ==2
    figure;
    visPolynomial2(coeffs, ws_bounds)

    % vizUnsafeSet([-1.2, -1, -0.7, 0.7])
    % vizUnsafeSet([0.5, 0.7, -0.7, 0.7])
    % vizUnsafeSet([-1, 0.5, 0.5, 0.7])
    % vizUnsafeSet([-1, 0.5, -0.7, -0.5])
    % vizUnsafeSet([-0.57, -0.53, -0.17, -0.13])
    % vizUnsafeSet([-0.57, -0.53, 0.28, 0.32])
    % vizInitSet([-0.8, -0.6, 0, 0.2], eta)
    
    % figure
    % comp_coeffs = [0.31400990334728096087, -0.32494080426818400253, 0.36104929238453525864, -0.07190328078594943051, -0.10862570595028053777, -0.03672242516433132931, -0.01768357616120186382, -0.01527958258179484829, 0.00240399357940715430]
    % visPolynomial2(comp_coeffs, ws_bounds)
end