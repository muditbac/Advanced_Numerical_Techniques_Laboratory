% This code solves a partial differential equation using BTCS (Backward 
% Time Central Space)
% 
% Author: Mudit Bachhawat
% Roll: 13MA20023
% Creation Date: 12th March, 2015
% Last Updated on: 12th March, 2015
% 
% 
% du/dt = k * d2u/dx2
% u(x,0) = sin pi*x
% u(0,t) = 0, u(1,t) = 0, t = 1

k=1;
dt = 0.1;
dx = 0.02;
r = k*dt/dx^2;

% TODO Create a variable convention.
u0 = 0;
un = 0;

x0 = 0;
t0 = 0;

xn = 1;
tn = 0.5;

x = x0:dx:xn;
t = t0:dt:tn;

j = length(x) - 1;
n = length(t) - 1;

u = zeros(j+1, n+1);

for i = 1:j,
    u(i,1) = sin(pi*x(i));
end

% Solution

for i = 2:n+1,
    ai_x = @(x)( r );
    bi_x = @(x)-(1+2*r);
    ci_x = @(x)( r );
    di_x = @(x)( -x );


    mat_A = zeros(j-1, j-1);
    vec_b = zeros(j-1,1);

    % Preparing trigonal matrix
    for g=1:(j),
        mat_A(g,g) = bi_x(u(g+1,i-1));
        mat_A(g,g+1) = ci_x(u(g+1,i-1));
        mat_A(g+1,g) = ai_x(u(g+1,i-1));

        vec_b(g) = di_x(u(g+1,i-1));
    end
    
    vec_b(1) = vec_b(1) - r*u0;
    vec_b(j-1) = vec_b(j-1) - r*un;

    mat_A = mat_A(1:j-1, 1:j-1);
    vec_b = vec_b(1:j-1);

%     disp(mat_A)

    u(2:j, i) = thomas_algorithm(mat_A, vec_b);
    
    
end


mesh(t,x,u)
xlabel('Time')
ylabel('X')
zlabel('U(x)')