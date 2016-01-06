function [ output_args ] = runge_kutta( x, Y, h, f_derivative )
%RUNGE_KUTTA Returns value of Y at x+h using classical runge kutta method.
%   Input:
%       x: x value
%       Y: Y value
%       h: step size
%       f_derivative: derivative function
%   Output:
%       Y1: output at x+h

    k_1 = f_derivative(x,Y);
    k_2 = f_derivative(x+0.5*h,Y+0.5*h*k_1);
    k_3 = f_derivative((x+0.5*h),(Y+0.5*h*k_2));
    k_4 = f_derivative((x+h),(Y+k_3*h));

    output_args = Y + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;  % main equation
end