function [x, fval, iterates] = simplex(c, A, b, x0, outp)
%SIMPLEX Linear programming with the simplex method.
%   [X,FVAL,ITERATES] = simplex(c,A,b,x0,outp) attempts
%   to solve the linear programming problem:
%        
%            min c'*x    subject to:   A*x = b, x >= 0 
%             x
%   
%   The starting point is set by x0, and this point has to be feasible. If
%   outp is passed as 'report', the function prints a detailed report at
%   every iteration. If report is passed as [], no report is printed. At
%   termination, X contains the optimal point, FVAL the optimal objective
%   function, and ITERATES is a matrix where each column i contains x at 
%   iteration i.
%   
%   NOTE: This function is written for educational purposes only. The
%         algorithm is neither quick nor robust, and should only be used on
%         very small problems. The implementation is meant to illustrate
%         Procedure 13.1 in Nocedal and Wright: Numerical Optimization,
%         2nd ed., Springer, 2006. There are probably some errors in the
%         implementation, please report any to joakim.r.andersen@ntnu.no.
% 
%   Version 0.9
%   Tor Aksel N. Heirung, 31.01.2011
%
%   Changed contact info
%   Andreas L. Fl?ten, 22.01.2018
%
%   Changed contact info
%   Joakim R. Andersen, 15.01.2019
    

n = numel(c);
m = numel(b);

I = eye(n); % Identity matrix

n_corners = nchoosek(n,m);
iterates = NaN(n,n_corners); % A matrix containing x for all iterations

basic = find(x0)'; % Indices of nonzero elements in x0
nbasic = find(x0==0)'; % Indices of zero elements in x0

if (numel(basic) ~= m) || (numel(nbasic) ~= n-m)
    error('Invalid x0. x0 must contain m nonzero elements.');
elseif A*x0 ~= b
    error('x0 is infeasible.')
end

k = 1; % Iteration counter
opt = false;
while 1
    % Init matrices
    B = A(:,basic);
    N = A(:,nbasic);
    xB = B\b;
    cB = c(basic);
    cN = c(nbasic);
    la = B'\cB;
    sN = cN - N'*la;
    
    P = I(:, [basic' nbasic']); % Creating permutation matrix
    x = P*[xB; zeros(n-m,1)];
    fval = c'*x;
    iterates(:,k) = x;
        
    if prod((sN >= 0) + 0) % adding 0 for logical -> numeric
        opt = true;
        iterates = iterates(:, 1:k);
        if strcmpi(outp, 'report')
            report(k, basic, nbasic, xB, la, sN, x, fval, q, d, x_q_plus, p, pB, xB_plus, xN_plus, opt); 
            finalreport(x, fval);
        end        
        break;
    end
    i_min_sN = find(sN == min(sN), 1, 'first');
    q = nbasic(i_min_sN);
    d = B\A(:,q);
    if prod((d < 0) + 0) % adding 0 for logical -> numeric
        error('Problem unbounded.');
    end
    i_d_g_0 = find(d>0); % Indices of the elemenents in d greater than 0
    ratios = xB(i_d_g_0)./d(i_d_g_0);
    x_q_plus = min(ratios);
    p = i_d_g_0(find(ratios == x_q_plus, 1, 'first'));
    pB = basic(p);
    xB_plus = xB - d*x_q_plus;
    xN_plus = zeros(n-m,1); xN_plus(nbasic(q)) = x_q_plus;
    
    if strcmpi(outp, 'report')
        report(k, basic, nbasic, xB, la, sN, x, fval, q, d, x_q_plus, p, pB, xB_plus, xN_plus, opt); 
    end
    
    [basic, nbasic] = newBasis(basic, nbasic, q, pB);
    k = k+1;
    if k >= n_corners
        error('Iteration limit reached. Cycling?');
    end
end

end

function [basic, nbasic] = newBasis(basic, nbasic, entering, leaving)
    iLeaving = find(basic == leaving, 1, 'first');
    iEnter = find(nbasic == entering, 1, 'first');
    basic(iLeaving) = entering;
    nbasic(iEnter) = leaving;
end

function [] = report(iteration, basic, nbasic, xB, lambda, sN, x, fval, q, d, x_q_plus, p, pB, xB_plus, xN_plus, opt)
    line = repmat('-',1,80);
    fprintf('\n%s\n',line);
    fprintf('    Iteration number: %i\n', iteration);
    fprintf(['    Basic index set:    {' sprintf('%i, ', basic) '\b\b}\n' ]);
    fprintf(['    Nonbasic index set: {' sprintf('%i, ', nbasic) '\b\b}\n' ]);
    fprintf(['    x_B =    [' sprintf('%9.4f, ', xB) '\b\b]''\n' ]);
    fprintf(['    x_N =    [' sprintf('%9.4f, ', zeros(numel(nbasic),1)) '\b\b]''\n' ]);
    fprintf(['    lambda = [' sprintf('%9.4f, ', lambda) '\b\b]''\n' ]);
    fprintf(['    s_N =    [' sprintf('%9.4f, ', sN) '\b\b]''\n' ]);
    if opt
        return;
    end
    fprintf(['    x =      [' sprintf('%9.4f, ', x) '\b\b]''\n' ]);
    fprintf( '    c''x = %14.8f\n', fval);
    fprintf( '    x_%i will enter the basis (q = %i)\n', q, q);
    fprintf(['    d =      [' sprintf('%9.4f, ', d) '\b\b]''\n' ]);
    fprintf( '    x_q+ = x_%i+ = %9.4f \t(value of entering variable/step length)\n', q, x_q_plus);
    fprintf( '    x_%i will leave the basis (p = %i)\n', pB, p);
    fprintf(['    x_B+ =   [' sprintf('%9.4f, ', xB_plus) '\b\b]'' \t(Current basic vector at new point)\n' ]);
    fprintf(['    x_N+ =   [' sprintf('%9.4f, ', xN_plus) '\b\b]'' \t(Current nonbasic vector at new point)\n' ]);
    fprintf('%s\n\n',line);
end

function [] = finalreport(x, fval)
    line = repmat('-',1,80);
    fprintf('\n');
    fprintf('    OPTIMAL POINT FOUND\n');
    fprintf(['    x^* =    [' sprintf('%9.4f, ', x) '\b\b]''\n' ]);
    fprintf( '    c''x^* = %14.8f\n', fval);
    fprintf('%s\n\n',line);
end