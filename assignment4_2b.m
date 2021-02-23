% Problem data: (Remember that all vectors are column vectors.)
c = [-3/2 -1]';
A = [-2   -1;
     -1 -3];
b = [-8 -15]';
% A fesible starting point:
x = linspace(-10, 10, 100);
y = linspace(-10, 10, 100);
[X, Y] = meshgrid(x, y);
Z = -(3-0.4*X).*X - (2-0.2*Y).*Y; 
levels = (-12:2:8)';
contour(X, Y, Z, levels, 'Color', .7*[1 1 1]);
hold on
c1 = 8 - 2*x;
c2 = (15 - x)/3;
plot(x, c1);
plot(x, c2);
x = [0 2.4 2.25];
y = [0 3.2 3.5];
labels = {'Iteration 0','Iteration 1','Iteration 2'};
plot(x,y,'o')
text(x,y,labels,'VerticalAlignment','bottom','HorizontalAlignment','right')
legend('Objective function f(x)','Constraint RI', 'Constraint RII');
xlim([-2, 10]);
ylim([-2, 10]);
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
hold off