clear; close all; clc;

N = 5;
d = 2;
ws_bounds = [0, 1, 0, 1];

file_directory = "../build/bin/";

A = readMatrixFromFile(file_directory + "A.txt");
b = readMatrixFromFile(file_directory + "b.txt");
Phi_inv = readMatrixFromFile(file_directory + "Phi_inv.txt");

objective_vec = zeros(size(A,2), 1);
objective_vec(end - 1) = 1.0;
objective_vec(end) = N;

vars = linprog(objective_vec, -A, -b);

beta = vars(1:end-2);
eta = vars(end - 1);

coeffs = Phi_inv * beta;
% disp(coeffs)

if d == 1
    visPolynomial(coeffs, "red")
elseif d ==2
    visPolynomial2(coeffs, ws_bounds)
end