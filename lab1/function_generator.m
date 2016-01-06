function [ r_fn ] = function_generator( F_xy )
%FUNCTION_GENERATOR Summary of this function goes here
%   Detailed explanation goes here




r_fn = @f_derivative;

    function [ output_args ] = f_derivative( x, Y )
    %DERIVATIVE Returns derivatives for nth order ODE
    %   Input
    %       x: Input point
    %       Y: Y values in order [y;y';y''.. ] till y(n-1)

    len = length(Y);

    F_Y = zeros(size(Y));

        for i = 1:len-1 
            F_Y(i) = Y(i+1);
        end

    F_Y(len) = F_xy(x, Y);

    output_args = F_Y;

    end



end

