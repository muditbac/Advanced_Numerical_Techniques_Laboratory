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
% mat_A -  4D tridiagonal matrix containing the coefficients of equations
% b - 3D RHS of the equations
% Returns:
% output_args - Solution of the equations

function [ output_args ] = thomas_algorithm_block( mat_A,b )
    
    l = length(mat_A);
    
    if l==2,
        output_args = mat_A\b;
        return
    end
    
    for i=1:l-1,
        fac = mat_A(:,:,i+1,i) / mat_A(:,:,i,i);
       
        mat_A(:,:,i+1,i) = mat_A(:,:,i+1,i) - fac * mat_A(:,:,i,i);
        mat_A(:,:,i+1,i+1) = mat_A(:,:,i+1,i+1) - fac * mat_A(:,:,i,i+1);
        
        if i~=l-1,
            mat_A(:,:,i+1,i+2) = mat_A(:,:,i+1,i+2) - fac * mat_A(:,:,i,i+2);
        end
        b(:,:,i+1,:) = b(:,:,i+1,:) - fac*b(:,:,i,:);
    end
    for i=l:-1:2,
        fac = mat_A(:,:,i-1,i) / mat_A(:,:,i,i);
        
        mat_A(:,:,i-1,i) = mat_A(:,:,i-1,i)- fac*mat_A(:,:,i,i);
        
        b(:,:,i-1,:) = b(:,:,i-1,:) - fac*b(:,:,i,:);
        
        output_args(:,i) = mat_A(:,:,i,i) \ b(:,:,i,:);
    end
    
    output_args(:,1) = mat_A(:,:,1,1) \ b(:,:,1,:);
    
%     output_args = b ./ diag(mat_A);

end

