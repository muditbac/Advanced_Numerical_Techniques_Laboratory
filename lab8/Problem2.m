% Solve Possions Equation using 5-Point Method and Gauss Seidal Method
% 
% 
% Author: Mudit Bachhawat
% Creation Date: 3rd April, 2016
% Last Updated on: 3rd April, 2016
% 
% 
% d2u/dx2 + d2u/dy2 = x^2 + y^2
% u(x ,y) = 0 for 0 <= x,y <= 1
% h = 1/4

clear, clc;
h = 0.25;
dx = h;
dy = h;

x = 0:dx:1;
y = 0:dy:1;

m = length(x)-1;
n = length(y)-1;

u(:,:,:) = zeros(m+1, n+1);
z = 1;
first = true;
while first || max(max(abs(u(:,:,end)-u(:,:,end-1))))>1e-5,
    first  = false;
    z = z+1;
    u(:,:,z) = zeros(m+1, n+1);
    for i = 2:m,
        for j = 2:n,
            u(i,j,z) = 1/4*(u(i+1,j,z-1)+u(i-1,j,z)+u(i,j+1,z-1)+u(i,j-1,z) - h*h*x(i)*x(i) - h*h*y(j)*y(j));
        end
    end
end

mesh(x,y,u(:,:,end));
xlabel('X');
ylabel('Y')
figure;

% Animation Gauss Seidal
for i = 1:length(u),
    mesh(x,y,u(:,:,i));
    axis([0 1 0 1 -0.05 0])
    drawnow
end

