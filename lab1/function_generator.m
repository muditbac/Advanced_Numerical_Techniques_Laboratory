% This code is used to convert derivative of nth order ODE to n, 1st order 
% ODEs. It inputs explicit equation of nth order derivative of equation.
% 
% Author: Mudit Bachhawat
% Creation Date: 5th January, 2015
% Last Updated on: 6th January, 2015
% 
% The meaning if input and output parameters:
% 
% Input:
% F_xy -  explicit equation of nth order derivative of equation.
% Returns:
% r_fn - returns function which inputs x and Y and returns derivative of
% y and all other derivatives.


function [ r_fn ] = function_generator( F_xy )

r_fn = @f_derivative;

    function [ output_args ] = f_derivative( x, Y )

    len = length(Y);

    F_Y = zeros(size(Y));

        for i = 1:len-1 
            % Pushes array and update last derivative
            F_Y(i) = Y(i+1);
        end

    F_Y(len) = F_xy(x, Y);

    output_args = F_Y;

    end

end

