clear; close all;

file_directory = "../../build/bin/";

p = readMatrixFromFile(file_directory + "certificate_coeffs.txt");
ws_bounds = [-1.2, .7, -.7, .7];
% ws_bounds = [0.0, .9, 0.0, .9];
% ws_bounds = [0, 1, 0, 1];

visPolynomial2(p, ws_bounds)
axis([ws_bounds, 0, 1])


eta = 0.03;

% vizUnsafeSet([-.57, -.53, -.17, -.13])
% vizUnsafeSet([-1.2, -1, -0.7, 0.7])
% vizUnsafeSet([0.5, 0.7, -0.7, 0.7])
% vizUnsafeSet([-1, 0.5, 0.5, 0.7])
% vizUnsafeSet([-1, 0.5, -0.7, -0.5])
vizUnsafeSet([-0.57, -0.53, -0.17, -0.13])
vizUnsafeSet([-0.57, -0.53, 0.28, 0.32])
vizInitSet([-0.8, -0.6, 0, 0.2], eta)




