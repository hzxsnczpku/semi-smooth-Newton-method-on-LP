function [x, out] = lp_drs_newton(c, A, b, opts, x0)
%  --------------------------------------------------------------
%  LP DRS Newton Method
%
%  This function solves the LP problem
%
%     minimize    c^Tx
%     subject to  Ax = b
%                 x >= 0
%
%  using the DRS Newton method. 
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
        t = 5;              
    end
    
    if isfield(opts, 'gamma1')
        gamma1 = opts.gamma1;
    else
        gamma1 = 1.1;
    end
    
    if isfield(opts, 'gamma2')
        gamma2 = opts.gamma2;
    else
        gamma2 = 1.5;
    end
    
    if isfield(opts, 'eta1')
        eta1 = opts.eta1;
    else
        eta1 = 0.1;
    end
    
    if isfield(opts, 'eta2')
        eta2 = opts.eta2;
    else
        eta2 = 0.8;
    end
    
    if isfield(opts, 'v')
        v = opts.v;
    else
        v = 0.4;
    end
    
    if isfield(opts, 'iters')   % outer loop
        iters = opts.iters;
    else
        iters = 3000;
    end
    
    if isfield(opts, 'err1')
        err1 = opts.err1;
    else
        err1 = 1e-8;
    end

    if isfield(opts, 'lam_')
        lam_ = opts.lam_;
    else
        lam_ = 1e-6;
    end
    
    if isfield(opts, 'lam')
        lam = opts.lam;
    else
        lam = 0.05;
    end

    %% Initialization
    [~, n] = size(A);

    P = A' * ((A * A') \ A);
    D = eye(n) - P;
    T = 2 * P - eye(n) ;
    B = A' * ( (A * A') \ b);

    F = @(z) T * max(z - t * c, 0) + D * z - B;
    dF = @(z) T * diag(z - t * c > 0) + D;

    k = 0;
    z = zeros(n, 1);
    u_ = z;
    fu_ = F(u_);
    out.phistory = [];

    %% Main Loop
    while k <= iters
        k = k + 1;
        fz = F(z);
        d = -(dF(z) + lam * norm(fz) * eye(n)) \ fz;

        dnorm = norm(d);
        if dnorm < err1
            break;
        end

        u = z + d;
        fu = F(u);
        dp = -fu' * d;
        rho = dp / dnorm ^ 2;

        if rho >= eta1
            if norm(fu) <= v * norm(fu_)
                z = u; 
                fu_ = fu;
            else
                z = z - dp / norm(fu) ^ 2 * fu;
            end
        end

        if rho >= eta2
            lam = (lam_ + lam) / 2;
        elseif rho >= eta1
            lam = (lam + gamma1 * lam) / 2;
        else
            lam = (gamma1 + gamma2) / 2 * lam;
        end

        out.phistory = [out.phistory, c' * max(z - t * c, 0)];
    end

    %% Output
    x = max(z - t * c, 0); 
    out.pobjval = c'* x;
    out.k = k;
    out.x = x;
    out.pfeasibility = norm(A * x - b);
    
end
