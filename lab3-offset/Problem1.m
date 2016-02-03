% This code solves a boundary value problem using finite difference and
% thomas algorithm
% 
% Author: Mudit Bachhawat
% Roll: 13MA20023
% Creation Date: 2nd Feburary, 2015
% Last Updated on: 2nd Feburary, 2015
% 
% Equation 
% y'''' + 81y = 81*x^2
% y0, yn, y''0, y''n are given

h = 0.05;

x0 = 0;
xn = 1;

y0 = 0;
yn = 0;

z0 = 0;
zn = 0;

n = (xn - x0)/h - 1;

x = x0+h:h:xn-h;
y(n) = 0;
z(n) = 0;
w = [y(1);z(1)];

for i=2:n
    w(:,:,i) = [y(i);z(i)];
end

A = [1/h.^2 0;0 1/h.^2];
for i=2:n
    A(:,:,i)= [1/h.^2 0;0 1/h.^2]; 
end

B = [-2/h.^2 -1;81 -2/h.^2];
for i=2:n
    B(:,:,i)= [-2/h.^2 -1;81 -2/h.^2]; 
end

C = [1/h.^2 0;0 1/h.^2];
for i=2:n
    C(:,:,i)= [1/h.^2 0;0 1/h.^2]; 
end

D = [0-A(1,1,1)*y0-A(1,2,1)*z0 ; 81*x(1).^2-A(2,1,1)*y0-A(2,2,1)*z0];

for i=2:n-1
    D(:,:,i)= [0;81*x(i).^2];
end

D(:,:,n) = [0-C(1,1,n)*yn-C(1,2,n)*zn ; 81*x(n).^2-C(2,1,n)*yn-C(2,2,n)*zn];

f_Cp(2,2) = 0;

% Solving TriDiagonal Matrix
for i=1:n-1
   if i==1
       f_Cp(:,:,i) = inv(B(:,:,i))*C(:,:,i);
   else
       f_Cp(:,:,i) = inv((B(:,:,i)-A(:,:,i)*f_Cp(:,:,i-1)))*C(:,:,i);
   end
end

f_Dp(2,1) = 0;

for i=1:n
   if i==1
       f_Dp(:,:,i) = inv(B(:,:,i))*D(:,:,i);
   else
       f_Dp(:,:,i) = inv((B(:,:,i)-A(:,:,i)*f_Cp(:,:,i-1)))*(D(:,:,i)-A(:,:,i)*f_Dp(:,:,i-1));
   end
end


for i=n:-1:1
    if i==n
        w(:,:,i) = f_Dp(:,:,i);
    else
        w(:,:,i) = f_Dp(:,:,i)-f_Cp(:,:,i)*w(:,:,i+1);
    end
end

for i=1:n
   y(i) = w(1,1,i);
end


plot(x,y);
xlabel('X')
ylabel('Y')
