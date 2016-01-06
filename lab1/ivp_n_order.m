% This code implements initial value problem for n, 1st order equations with
% same x variable.
% 
% Author: Mudit Bachhawat
% Creation Date: 5th January, 2015
% Last Updated on: 6th January, 2015
% 
% The meaning if input and output parameters:
% 
% Input:
% x_val - inital point
% Y_val - values of all variables(if any), 
% h - step size
% fn_derivative - function taking input as x and Y in specified format,
%                 returns derivatives of all unknown variables
% max_x - maximum x value to which solution has to be calulated
% 
% Returns:
% x - x values whose y values are calculated
% y - values of all the variables.
% 
% The example of using:
% fexp = @(x,y)(x);
% [x,y] = ivp_n_order(0,0,0.1, fexp , 1)
% It is meant that will be solved the IVP ODE y' = x described in the function fun,
% on the interval (0,1) with boundary conditions y(0) = 0


function [ x, y ] = ivp_n_order( x_val, Y_val, h, fn_derivative, max_x )

x = x_val:h:max_x;
l = length(Y_val);
y = zeros(l,length(x)); 
y(:,1) = Y_val; % initial condition

for i=1:(length(x)-1)  
%     Classical Runge Kutta Method
    k_1 = fn_derivative(x(i),y(:,i));
    k_2 = fn_derivative(x(i)+0.5*h,y(:,i)+0.5*h*k_1);
    k_3 = fn_derivative((x(i)+0.5*h),(y(:,i)+0.5*h*k_2));
    k_4 = fn_derivative((x(i)+h),(y(:,i)+k_3*h));
    
    y(:,i+1) = y(:,i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;
end


end

