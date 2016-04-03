% This function solves a tridiagonal set of equations using thomas
% algorithm.
% 
% Author: Mudit Bachhawat
% Creation Date: 19th January, 2015
% Last Updated on: 19th January, 2015
% 
% The meaning if input and output parameters:
% 
% Input:
% mat_A -  tridiagonal matrix containing the coefficients of equations
% b - RHS of the equations
% Returns:
% output_args - Solution of the equations

function [ output_args ] = thomas_algorithm( mat_A,b )
    
    l = length(mat_A);
    
    for i=1:l-1,
        fac = mat_A(i+1,i) / mat_A(i,i);
        mat_A(i+1,:) = mat_A(i+1,:) - fac*mat_A(i,:);
        b(i+1,:) = b(i+1,:) - fac*b(i,:);
    end
    for i=l:-1:2,
        fac = mat_A(i-1,i) / mat_A(i,i);
        mat_A(i-1,:) = mat_A(i-1,:) - fac*mat_A(i,:);
        b(i-1,:) = b(i-1,:) - fac*b(i,:);
    end
    
    output_args = b ./ diag(mat_A);

end

