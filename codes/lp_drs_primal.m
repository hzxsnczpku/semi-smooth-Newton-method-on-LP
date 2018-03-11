function [x, out] = lp_drs_primal(c, A, b, opts, x0)
%  --------------------------------------------------------------
%  LP DRS Method
%
%  This function solves the LP problem
%
%     minimize    c^Tx
%     subject to  Ax = b
%                 x >= 0
%
%  using the DRS method. 
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
    if isfield(opts, 't')   % parameter for the penalty term          
        t = opts.t;
    else
        t = 10;              
    end
    
    if isfield(opts, 'gamma')
        gamma = opts.gamma;
    else
        gamma = 1;
    end
    
    if isfield(opts, 'iters')   % outer loop
        iters = opts.iters;
    else
        iters = 20000;
    end
    
    if isfield(opts, 'err1')
        err1 = opts.err1;
    else
        err1 = 1e-8;
    end
    
    if isfield(opts, 'err2')
        err2 = opts.err1;
    else
        err2 = 1e-8;
    end

    %% Initialization
    [~, n] = size(A);
    
    k = 0;
    w = randn(n, 1);
    x = x0;
    D = eye(n) - A' * ((A * A') \ A);
    B = A' * ((A*A') \ b);
    
    out.phistory = [];

    %% Main Loop
    while k < iters
        k = k + 1;
        u = max(x - w - t * c, 0);   
        x = D * (u + w) + B;
        w = w + gamma * (u - x);

        if mod(k, 100) == 0
            if norm(u - x) < err1 && norm(A * max(x, 0) - b) < err2
                break;
            end
        end

        out.phistory = [out.phistory, c' * max(x, 0)];
    end

    %% Output
    x = max(x, 0); 
    out.pobjval = c' * x;
    out.k = k;
    out.x = x;
    out.pfeasibility = norm(A * x - b);

end
