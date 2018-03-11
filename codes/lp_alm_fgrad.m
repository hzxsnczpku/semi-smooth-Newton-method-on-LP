function [x, out] = lp_alm_fgrad(c, A, b, opts, x0)
%  --------------------------------------------------------------
%  LP ALM Fast Gradient Method
%
%  This function solves the LP problem
%
%     minimize    c^Tx
%     subject to  Ax = b
%                 x >= 0
%
%  using the ALM method with the fast gradient method. 
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
        iters = 500;
    end
    
    if isfield(opts, 'inner_iters')   % inner loop
        inner_iters = opts.inner_iters;
    else
        inner_iters = 100;
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
    
    if isfield(opts, 'err')
        err = opts.err;
    else
        err = 1e-8;
    end
    
    if isfield(opts, 'alpha')    % the step size
        alpha = opts.alpha;
    else
        alpha = 8e-4;
    end

    %% Initialization
    [m, ~] = size(A);
    
    k = 0;
    y = randn(m,1);
    x = x0 / t;
    
    out.phistory = [];
    out.dhistory = [];

    %% Main Loop
    while k < iters
        k = k + 1;
        y0 = y;
        for i = 1 : inner_iters
            z = y + (i - 2) / (i + 1) * (y - y0); 
            g = -b + t * (A * max(A' * z + x - c, 0));
            if norm(g) < err
                break
            end
            y0 = y;
            y = z - alpha * g;
        end

        Ay = A' * y;
        s = max(c - x - Ay, 0);
        dx = gamma * (Ay + s - c);
        x = x + dx;
        x_t = x * t;

        out.phistory = [out.phistory, c' * max(x_t, 0)];
        out.dhistory = [out.dhistory, b'* y];

        if norm(dx) < err1 && norm(A * x_t - b) < err2 
            break
        end
    end

    %% Output
    x = max(x_t, 0); 
    out.pobjval = c'* x;
    out.dobjval = b'* y;
    out.k = k;
    out.x = x;
    out.s = s;
    out.y = y;
    out.pfeasibility = norm(A * x - b);
    out.dfeasibility = norm(A' * y + s - c);
end
