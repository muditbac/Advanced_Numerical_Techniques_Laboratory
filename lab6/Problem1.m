% This code solves a partial differential equation using Crank Nicolson
% Method
% 
% Author: Mudit Bachhawat
% Roll: 13MA20023
% Creation Date: 12th March, 2016
% Last Updated on: 12th March, 2016
% 
% 
% du/dt = k * d2u/dx2
% u(x,0) = cos pi*x/2
% u(-1,t) = 0, u(1,t) = 0, t = 1

k=1;
dt = 1/27;
dx = 1/3;
r = k*dt/dx^2;

u0 = 0;
un = 0;

x0 = -1;
t0 = 0;

xn = 1;
tn = 0.8;

x = x0:dx:xn;
t = t0:dt:tn;

j = length(x) - 1;
n = length(t) - 1;

u = zeros(j+1, n+1);

for i = 1:j,
    u(i,1) = cos(pi*x(i)/2);
end

% Solution

ai_x = @(x)( -k/(2*dx*dx) );
bi_x = @(x)(1/dt + k/(dx*dx));
ci_x = @(x)( -k/(2*dx*dx) );
di_x = @(x,um1,u,up1)( (up1-2*u+um1)/(dx*dx) + u/dt );

for i = 2:n+1,

    mat_A = zeros(j-1, j-1);
    vec_b = zeros(j-1,1);

    % Preparing trigonal matrix
    for g=1:(j),
        mat_A(g,g) = bi_x(NaN);
        mat_A(g,g+1) = ci_x(NaN);
        mat_A(g+1,g) = ai_x(NaN);

        if (g~=j),
            vec_b(g) = di_x(NaN, u(g,i-1), u(g+1,i-1), u(g+2,i-1));
        end
    end
    
    
    
    vec_b(1) = vec_b(1) - r*u0;
    vec_b(j-1) = vec_b(j-1) - r*un;

    mat_A = mat_A(1:j-1, 1:j-1);

    vec_b = vec_b(1:j-1);

    u(2:j, i) = thomas_algorithm(mat_A, vec_b);
    
    
end


mesh(t,x,u)
xlabel('Time')
ylabel('X')
zlabel('U(x)')