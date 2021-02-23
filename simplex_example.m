% Example illustrating the use of the function simplex. We here solve
% Example 13.1 in the textbook. (Note that this example has many errors in 
% the text.)

clear('all');
clc;

% Problem data: (Remember that all vectors are column vectors.)
c = [-4 -2 0 0]';
A = [1   1  1 0 ;
     2 1/2  0 1];
b = [5 8]';
% A fesible starting point:
x0 = [0 0 5 8]'; % The starting point must also have four variables!

% Solve problem and print a report: (Call simplex(c,A,b,x0,[]) if you do 
% not want the report printed).
[x, fval, iterates] = simplex(c,A,b,x0,'report');

% Extract iterates in the space of x_1 and x_2:
iter_x1_x2 = iterates(1:2, :);

% Extract iterates as individual vectors (containing x_1 and x_2):
iter_1 = iter_x1_x2(:,1);
iter_2 = iter_x1_x2(:,2);
iter_3 = iter_x1_x2(:,3);