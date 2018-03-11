% function Test_lp_dual

%% random data
rng(5)
% rng(2460)
% rng(20)

n = 100;
m = 20;
A = sparse(rand(m, n));
xs = full(abs(sprandn(n, 1, m / n)));
b = A * xs;
y = randn(m, 1);
s = rand(n, 1) .* (xs == 0);
c = A'*y + s;
x0 = abs(randn(n, 1));

%% netlib
% load share2bpre.mat;
% out = preprocess(Model);
% c = out.c;
% A = out.A;
% b = out.b;
% [m, n] = size(A);
% x0 = abs(randn(n, 1));

%% calling cvx_mosek
opts1 = [];
tic;
[x1, out1] = lp_cvx_mosek(c, A, b, opts1, x0);
t1 = toc;

%% calling cvx_gurobi
opts2 = [];
tic;
[x2, out2] = lp_cvx_gurobi(c, A, b, opts2, x0);
t2 = toc;

%% calling mosek
opts3 = [];
tic;
[x3, out3] = lp_mosek(c, A, b, opts3, x0);
t3 = toc;

%% calling gurobi
opts4 = [];
tic;
[x4, out4] = lp_gurobi(c, A, b, opts4, x0);
t4 = toc;

%% calling ALM
opts5 = []; 
tic;
[x5, out5] = lp_alm(c, A, b, opts5, x0);
t5 = toc;

%% calling fast gradient ALM 
opts6 = []; 
tic;
[x6, out6] = lp_alm_fgrad(c, A, b, opts6, x0);
t6 = toc;

%% calling semi-smooth newton ALM
opts7 = []; 
tic;
[x7, out7] = lp_alm_newton(c, A, b, opts7, x0);
t7 = toc;

%% calling ADMM
opts8 = []; 
tic;
[x8, out8] = lp_admm_dual(c, A, b, opts8, x0);
t8 = toc;

%% calling DRS
opts9 = []; 
tic;
[x9, out9] = lp_drs_primal(c, A, b, opts9, x0);
t9 = toc;
 
%% calling DRS Newton
opts10 = []; 
tic;
[x10, out10] = lp_drs_newton(c, A, b, opts10, x0);
t10 = toc;

%% print comparison results with mosek lp solvers
errfun = @(x1, x2) norm(x1(c~=0)-x2(c~=0))/(1 + norm(x1(c~=0)));
fprintf('cvx-call-mosek:      obj: %5.12f, cpu: %5.2f, err-to-cvx-mosek: %3.2e\n', out1.objval,  t1, errfun(x1, x1));
fprintf('cvx-call-gurobi:     obj: %5.12f, cpu: %5.2f, err-to-cvx-mosek: %3.2e\n', out2.objval,  t2, errfun(x1, x2));
fprintf('call-mosek:          obj: %5.12f, cpu: %5.2f, err-to-cvx-mosek: %3.2e\n', out3.pobjval, t3, errfun(x1, x3));
fprintf('call-gurobi:         obj: %5.12f, cpu: %5.2f, err-to-cvx-mosek: %3.2e\n', out4.objval,  t4, errfun(x1, x4));
fprintf('call-ALM:            obj: %5.12f, cpu: %5.2f, err-to-cvx-mosek: %3.2e\n', out5.pobjval, t5, errfun(x1, x5));
fprintf('call-ALM_fgrad:      obj: %5.12f, cpu: %5.2f, err-to-cvx-mosek: %3.2e\n', out6.pobjval, t6, errfun(x1, x6));
fprintf('call-ALM_Newton:     obj: %5.12f, cpu: %5.2f, err-to-cvx-mosek: %3.2e\n', out7.pobjval, t7, errfun(x1, x7));
fprintf('call-ADMM:           obj: %5.12f, cpu: %5.2f, err-to-cvx-mosek: %3.2e\n', out8.pobjval, t8, errfun(x1, x8));
fprintf('call-DRS:            obj: %5.12f, cpu: %5.2f, err-to-cvx-mosek: %3.2e\n', out9.pobjval, t9, errfun(x1, x9));
fprintf('call-DRS_Newton:     obj: %5.12f, cpu: %5.2f, err-to-cvx-mosek: %3.2e\n', out10.pobjval, t10, errfun(x1, x10));

% subplot(5,2,1);plot(x1);
% subplot(5,2,2);plot(x2);
% subplot(5,2,3);plot(x3);
% subplot(5,2,4);plot(x4);
% subplot(5,2,5);plot(x5);
% subplot(5,2,6);plot(x6);
% subplot(5,2,7);plot(x7);
% subplot(5,2,8);plot(x8);
% subplot(5,2,9);plot(x9);
% subplot(5,2,10);plot(x10);
