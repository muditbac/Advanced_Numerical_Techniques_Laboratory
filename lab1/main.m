h=0.01;                                             % step size
x = 0:h:6;                                         % Calculates upto y(3)
y = zeros(2,length(x)); 
y(:,1) = [1;0];                                          % initial condition


% Here that Y is array of all y values in order [y y' y'' ...]
% F_xy = @(x,Y) Y(1);                    % change the function as you desire
for i=1:(length(x)-1)                              % calculation loop
    
    
    k_1 = f_derivative(x(i),y(:,i));
    k_2 = f_derivative(x(i)+0.5*h,y(:,i)+0.5*h*k_1);
    k_3 = f_derivative((x(i)+0.5*h),(y(:,i)+0.5*h*k_2));
    k_4 = f_derivative((x(i)+h),(y(:,i)+k_3*h));
    
    y(:,i+1) = y(:,i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;
end

plot (x,y)