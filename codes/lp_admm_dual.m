function [x, out] = lp_admm_dual(c, A, b, opts, x0)
%  --------------------------------------------------------------
%  LP Dual ADMM Method
%
%  This function solves the dual LP problem
%
%     minimize    b^Ty
%     subject to  A^Ty + s = c
%                 s >= 0
%
%  using the the ADMM method. 
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

    %% Hyperparameters
    if isfield(opts, 't')              
        t = opts.t;
    else
        t = 20;              
    end
    
    if isfield(opts, 'gamma')              
        gamma = opts.gamma;
    else
        gamma = 1.618;              
    end
    
    if isfield(opts, 'iters')              
        iters = opts.iters;
    else
        iters = 2e4;              
    end
    
    if isfield(opts, 'err1')              
        err1 = opts.err1;
    else
        err1 = 1e-10;              
    end
    
    if isfield(opts, 'err2')              
        err2 = opts.err2;
    else
        err2 = 1e-10;              
    end

    %% Initialization
    [m, n] = size(A);

    k = 0;
    y = randn(m, 1);
    x = x0 / t;
    s = zeros(n, 1);
    L = chol(t * (A * A'), 'lower');
    Ay = A' * y; 
    
    out.phistory = [];
    out.dhistory = [];

    %% Main Loop
    while k <= iters
        k = k + 1;
        s = max(c - x - Ay, 0);
        y = L' \ (L \ (t * (A * (c - s - x)) + b));
        Ay = A' * y; 
        dx = gamma * (Ay + s - c);
        x = x + dx;
        x = max(x, 0);
        x_t = t * x;

        if norm(dx) < err1 && norm(A * max(x_t, 0) -b) < err2
            break;
        end

        out.phistory = [out.phistory, c' * max(x_t, 0)];
        out.dhistory = [out.dhistory, b'* y];
    end

    %% Output
    x = max(x_t, 0); 
    out.pobjval = c'* x;
    out.dobjval = b'* y;
    out.k = k;
    out.x = x;
    out.y = y;
    out.s = s;
    out.pfeasibility = norm(A * x - b);
    out.dfeasibility = norm(A' * y + s - c);
    
end
