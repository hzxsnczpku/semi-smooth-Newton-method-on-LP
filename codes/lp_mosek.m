function [x, out] = lp_mosek(c, A, b, opts, x0)
%  --------------------------------------------------------------
%  LP Mosek
%
%  This function solves the LP problem
%
%     minimize    c^Tx
%     subject to  Ax = b
%                 x >= 0
%
%  calling Mosek directly. 
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
blc = b;
buc = b;
blx = zeros(size(x0));
bux = [];

%% Call Mosek
res = msklpopt(c, A, blc, buc, blx, bux);

%% Output
out = res.sol.itr;
x = out.xx;

end