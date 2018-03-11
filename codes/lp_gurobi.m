function [x, out] = lp_gurobi(c, A, b, opts, x0)
%  --------------------------------------------------------------
%  LP Gurobi
%
%  This function solves the LP problem
%
%     minimize    c^Tx
%     subject to  Ax = b
%                 x >= 0
%
%  calling Gurobi directly. 
%
%  Authors: Li Xiang,
%           Lin Dachao,
%           Ni Chengzhuo,
%           School of Mathematical Science, PKU
%  --------------------------------------------------------------
%
%  =========================== Inputs ===========================
%  
%     c: n * 1 matrix, the given vector of the object
%
%     A: m * n matrix, the given matrix of the constraint
%
%     b: m * 1 matrix, the given vector of the constraint
%
%  opts:    structure, modify options
%
%    x0: n * 1 matrix, the starting point of the algotirhm
%
%  ==============================================================
%
%  =========================== Outputs ==========================
%  
%     x: m * 1 matrix, the optimal point found by the algorithm
%
%   out:    structure, the record of the process information
%
%  ==============================================================

%% Initialization
[~, n] = size(A);
model.A = [A; -A; speye(n)];
model.obj = full(c)';
model.modelsense = 'Min';
model.rhs = full([b', -b', sparse(1, n)]);
model.sense = '>';

%% Call Gurobi
out = gurobi(model);

%% Output
x = out.x;

end

