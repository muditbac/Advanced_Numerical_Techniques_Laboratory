% This code solves a boundary value problem using finite difference method
% and thomas algorithm
% 
% Author: Mudit Bachhawat
% Roll: 13MA20023
% Creation Date: 12th January, 2015
% Last Updated on: 19th January, 2015
% 
% y'' + 2xy' + 2y = 4x
Ax = @(x)(2*x);
Bx = @(x)(2);
Cx = @(x)(4*x);

h = 0.1;

x0 = 0;
y0 = 1;

xn = 0.5;
yn = 1.279;

% Solution

ai_x = @(x)( 1/(h*h) - Ax(x)/(2*h) );
bi_x = @(x)( -2/(h*h) + Bx(x) );
ci_x = @(x)( 1/(h*h) + Ax(x)/(2*h) );
di_x = @(x)( Cx(x) );

x = x0:h:xn;

l = length(x);

mat_A = zeros(l-2,l-2);
vec_b = zeros(l-2,1);

% Preparing trigonal matrix
for i=1:(l-2),
    mat_A(i,i) = bi_x(x(i+1));
    mat_A(i,i+1) = ci_x(x(i+1));
    mat_A(i+1,i) = ai_x(x(i+2));
    
    vec_b(i) = di_x(x(i+1));
end

vec_b(1) = vec_b(1) - ai_x(x(1+1)) * y0;
vec_b(l-2) = vec_b(l-2) - ci_x(x((l-2) + 1)) * yn;
disp(mat_A)

mat_A = mat_A(1:l-2,1:l-2);
y = [y0;zeros(l-2,1);yn];

% Solves tridigonal equations using Thomas Algorithm
y(2:l-1) = thomas_algorithm(mat_A, vec_b);

disp(y(:));

plot(x,y)
legend('Solution for y = f(x)')
xlabel('X');
ylabel('Y');