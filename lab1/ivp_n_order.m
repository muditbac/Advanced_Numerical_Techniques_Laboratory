function [ x, y ] = ivp_n_order( x_val, Y_val, h, fn_derivative, max_x )
%IVP_STEP Summary of this function goes here
%   Documentation Remaining

x = x_val:h:max_x;                                         % Calculates upto y(3)
l = length(Y_val);
y = zeros(l,length(x)); 
y(:,1) = Y_val;                                          % initial condition

for i=1:(length(x)-1)                              % calculation loop
    k_1 = fn_derivative(x(i),y(:,i));
    k_2 = fn_derivative(x(i)+0.5*h,y(:,i)+0.5*h*k_1);
    k_3 = fn_derivative((x(i)+0.5*h),(y(:,i)+0.5*h*k_2));
    k_4 = fn_derivative((x(i)+h),(y(:,i)+k_3*h));
    
    y(:,i+1) = y(:,i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;
end

% plot (x,y)

end

