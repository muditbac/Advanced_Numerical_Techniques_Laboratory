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
% u(x,0) = 1, 0 < x < 1
% u(0,t) = u'(0,t), u(1,t) = -u'(1,t), t = 1

k=1;
dt = 0.04;
dx = 0.2;
r = k*dt/dx^2;

u0 = 1;
un = 1;

x0 = 0;
t0 = 0;

xn = 1;
tn = 1;

x = x0:dx:xn;
t = t0:dt:tn;

j = length(x) - 1;
n = length(t) - 1;

u = zeros(j+3, n+1);

for i = 1:j+3,
    u(i,1) = 1;
end

% Solution

ai_x = @(x)( -k/(2*dx*dx) );
bi_x = @(x)(1/dt + k/(dx*dx));
ci_x = @(x)( -k/(2*dx*dx) );
di_x = @(x,um1,u,up1)( (up1-2*u+um1)/(dx*dx) + u/dt );

for i = 2:n+1,

    mat_A = zeros(j+1, j+1);
    vec_b = zeros(j+1,1);

    % Preparing trigonal matrix
    for g=2:(j+2),
        g = g-1;
        mat_A(g,g) = bi_x(NaN);
        mat_A(g,g+1) = ci_x(NaN);
        mat_A(g+1,g) = ai_x(NaN);

        if (1),
            vec_b(g) = di_x(NaN, u(g,i-1), u(g+1,i-1), u(g+2,i-1));
        end
        g = g+1;
    end
    
    mat_A = mat_A(1:j+1, 1:j+1);
    
    mat_A(1,1) = mat_A(1,1) - 2*ai_x(NaN)*dx;
    mat_A(1,2) = mat_A(1,1) + ci_x(NaN);
    
    mat_A(end, end) = mat_A(end, end) - 2*dx*2*ci_x(NaN);
    mat_A(end, end-1) = mat_A(end, end-1) + ci_x(NaN);


    u(2:j+2, i) = thomas_algorithm(mat_A, vec_b);
    u(1, i) = u(3, i) - 2*dx*u(2,i);
    u(j+3, i) = u(end-2, i) - 2*dx*u(end-1,i);
%     u(1, i) = u0;
%     u(end, i) = un;
    
end

rotate3d on
mesh(t,x,u(2:end-1,:))
xlabel('Time')
ylabel('X')
zlabel('U(x)')