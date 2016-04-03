function [ output_args ] = spline_interpolation( x,y,y2, xx )
%SPLINE_INTERPOLATION Summary of this function goes here
%   Detailed explanation goes here
n = length(x);

yy = zeros(size(xx));

for i = 1:n-1,
    h = x(i+1)-x(i);
    f = @(xv)(y2(i)/6*((x(i+1)-xv)^3/h-h*(x(i+1)-xv)) + y(i)/h*(x(i+1)-xv) + y2(i+1)/6*((xv - x(i))^3/h-h*(xv - x(i))) + y(i+1)/h*(xv - x(i)));
    yy(xx<=x(i+1) & xx>x(i)) = arrayfun(f, xx(xx<=x(i+1) & xx>x(i)));
end
    
    output_args = yy;
end

