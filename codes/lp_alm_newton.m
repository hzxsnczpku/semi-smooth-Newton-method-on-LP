function [x, out] = lp_alm_newton(c, A, b, opts, x0)
%  --------------------------------------------------------------
%  LP ALM Newton Method
%
%  This function solves the LP problem
%
%     minimize    c^Tx
%     subject to  Ax = b
%                 x >= 0
%
%  using the ALM Newton method. 
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
    if isfield(opts, 'sig')   % parameter for the penalty term          
        sig = opts.sig;
    else
        sig = 1500;              
    end
    
    if isfield(opts, 'tau')
        tau = opts.tau;
    else
        tau = 1e-4;
    end
    
    if isfield(opts, 'tau1')
        tau1 = opts.tau1;
    else
        tau1 = 1e-4;
    end
    
    if isfield(opts, 'tau2')
        tau2 = opts.tau2;
    else
        tau2 = 1e-2;
    end
    
    if isfield(opts, 'iters')   % outer loop
        iters = opts.iters;
    else
        iters = 1000;
    end
    
    if isfield(opts, 'inner_iters')   % inner loop
        inner_iters = opts.inner_iters;
    else
        inner_iters = 50;
    end
    
    if isfield(opts, 'err')
        err = opts.err;
    else
        err = 1e-6;
    end

    if isfield(opts, 'mu')
        mu = opts.mu;
    else
        mu = 0.1;
    end
    
    if isfield(opts, 'delta')
        delta = opts.delta;
    else
        delta = 0.8;
    end

    %% Initialization
    c = -c;
    b = -b;
    A = -A;
    [m, ~] = size(A);

    phi = @(y, x) b' * y + 1 / 2 / sig * norm(max(x - sig * (A' * y - c), 0))^2;

    k = 0;
    x = x0;
    y = randn(m,1);
    out.phistory = [];
    out.dhistory = [];

    %% Main Loop
    while k <= iters
        k = k + 1;
        for i = 1:inner_iters 
            temp = max(x - sig * (A' * y - c), 0);
            g = b - A * temp;
            normg = norm(g);
            AA = A(:, (temp ~= 0));
            ep = tau1 * min(tau2, normg);
            h = sig * (AA * AA') + ep * speye(m);
            L = ichol(h, struct('type', 'ict', 'droptol', 1e-10, 'diagcomp', 0.0001));
            tol = min(1e-6, normg ^ (1 + tau));
            itr = 100;
            [d, ~] = pcg(h, g, tol, itr, L, L'); 
            alpha = 1;
            while phi(y - alpha * d, x) > phi(y, x) - mu * alpha * (g' * d)
                alpha = alpha * delta;
            end
            y = y - alpha * d;
            if normg < 1e-8
                break
            end
        end
        
        w = x - sig * (A' * y - c);
        x = max(w, 0);
        pobj = -c' * x;
        dobj = -b' * y;

        out.phistory = [out.phistory, pobj];
        out.dhistory = [out.dhistory, dobj];

        s = (x - w) / sig;
        RP = norm(A * x - b) / (1 + norm(b));
        RD = norm(-A' * y + s + c) / (1 + norm(c)); 
        gap = abs(pobj - dobj) /(1 + abs(pobj) + abs(dobj));
        if max([RP, RD, gap]) < err
            break
        end
        
    end

    %% Output
    out.pobjval = -c'* x;
    out.dobjval = -b'* y;
    out.k = k;
    out.x = x;
    out.y = y;
    out.s = s;
    out.pfeasibility = norm(A * x - b);
    out.dfeasibility = norm(A' * y - s - c);
    
end
