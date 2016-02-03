% This code solves a boundary value problem using shooting method
% 
% Author: Mudit Bachhawat
% Roll: 13MA20023
% Creation Date: 5th January, 2015
% Last Updated on: 6th January, 2015
% 
% Question 
% y' * y'' + 1 + (y')^2 = 0
% y(0) = 1
% y(1) = 2
% a0 = 0.5;
% a1 = 1;
% h = 0.2;
% expected result (alpha) = 2.0

x0 = 0;
xn = 1;

y0 = 1;
yn = 2;

a0 = 0.5;
a1 = 1;


h = 0.2;

d = @(x,Y) ((-Y(2)^2 -1)/Y(1));

% Generates derivative function for IVP
f_derivative = function_generator(d);

legend_vals = {};

while abs(a0-a1)>exp(-7)
    [~, y_vals]  = ivp_n_order(x0, [y0;a0], h, f_derivative, xn);

    disp(strcat('a0=', num2str(a0)))
    disp(y_vals)
    predicted_y_a0 = y_vals(1,end);

    [x_vals, y_vals]  = ivp_n_order(x0, [y0;a1], h, f_derivative, xn);

    disp(strcat('a1=', num2str(a1)))
    disp(y_vals)
    predicted_y_a1 = y_vals(1,end);

    a_new = a1 - (a1-a0)*(predicted_y_a1 - yn)/(predicted_y_a1 - predicted_y_a0);

    disp(a_new)

    a0 = a1;
    a1 = a_new;

    p = plot(x_vals, y_vals(1,:),'-.x', 'LineWidth',1, 'MarkerSize',5);
    legend_vals = [legend_vals (strcat('alpha = ',num2str(a0)))];
    hold on;
end

hold off;

set(p, 'LineWidth',1.5, 'LineStyle', '-')
xlabel('X Axis');
ylabel('Y Axis');

legend(legend_vals)