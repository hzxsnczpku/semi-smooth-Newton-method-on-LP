function [x, out] = lp_cvx_gurobi(c, A, b, opts, x0)
%  --------------------------------------------------------------
%  LP CVX Gurobi
%
%  This function solves the LP problem
%
%     minimize    c^Tx
%     subject to  Ax = b
%                 x >= 0
%
%  calling CVX Gurobi. 
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

    %% Problem Setting
    n = size(A, 2);

    cvx_begin quiet
    cvx_precision high
    cvx_solver gurobi
    variable x(n) 
    minimize(c' * x) 
    subject to
        A * x == b
        x >= 0
    cvx_end

    %% Output
    x = max(x, 0);
    out.status = 'OPTIMAL';
    out.objval = c' * x ;
    
end
