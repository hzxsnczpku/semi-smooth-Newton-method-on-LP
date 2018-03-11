function out = preprocess(s)
%% model standardization 
% from the problem                    
%   min   c'*x    
%   st.   bl <= A*x <= bu 
%         tl <=  x  <= tu
% equivalent to solve the problem  
%   min c'*x1-c'*x2
%   st. A(x1-x2)+s1 = bu
%       A(x1-x2)-s2 = bl
%         x1-x2 +s3 = tu
%         x1-x2 -s4 = tl
%       x1,x2,s1,s2,s3,s4 >= 0
[m,n] = size(s.A);
b = [s.rhs ; s.lhs ; s.ub ; s.lb];
A = [sparse(s.A), -sparse(s.A), speye(m),            sparse(m, m + 2 * n);
     sparse(s.A), -sparse(s.A), sparse(m, m), -speye(m), sparse(m, 2 * n);
     speye(n),    -speye(n),    sparse(n, 2 * m),  speye(n), sparse(n, n);
     speye(n),    -speye(n),    sparse(n, 2 * m), sparse(n, n), -speye(n)];
out.A = A(abs(b) ~= Inf, :);
out.b = b(abs(b) ~= Inf);
out.c = [s.obj; -s.obj; sparse(2 * m + 2 * n, 1)];
end

