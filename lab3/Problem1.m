% This code solves a boundary value problem using thomas algorithm
% 
% Author: Mudit Bachhawat
% Roll: 13MA20023
% Creation Date: 28th January, 2015
% Last Updated on: 28th January, 2015
% 
% Equation type
% y'' + Ay' + By = C
% a1y0 + b1y'0  = g1
% a2yn + b2y'n  = g2


Ax = @(x)(-2*x);
Bx = @(x)(-2);
Cx = @(x)(-4*x);

h = 0.1;

x0 = 0;
alpha1 = 1;
beta1 = -1;
gamma1 = 0;

xn = 1;
alpha2 = 2;
beta2 = -1;
gamma2 = 1;

% Solution

ai_x = @(x)( 1/(h*h) - Ax(x)/(2*h) );
bi_x = @(x)( -2/(h*h) + Bx(x) );
ci_x = @(x)( 1/(h*h) + Ax(x)/(2*h) );
di_x = @(x)( Cx(x) );

x = x0:h:xn;

l = length(x);

mat_A = zeros(l,l);
vec_b = zeros(l,1);

% Preparing trigonal matrix
for i=1:(l-1),
    mat_A(i,i) = bi_x(x(i));
    mat_A(i,i+1) = ci_x(x(i));
    mat_A(i+1,i) = ai_x(x(i+1));
    
    vec_b(i) = di_x(x(i));
end
% Compiling last element
mat_A(l, l) = bi_x(x(l));
vec_b(l) = di_x(x(l));

vec_b(1) = vec_b(1) + 2* ai_x(x(1)) * gamma1 * h / beta1;
vec_b(l) = vec_b(l) + 2* ci_x(x(l-1)) * gamma2 * h / beta2;


mat_A(1,1) = mat_A(1,1) + 2 * ai_x(x(1)) * alpha1 * h / beta1;
mat_A(1,2) = mat_A(1,2) + ci_x(x(1));

mat_A(l,l) = mat_A(l,l) + 2 * ci_x(x(l)) * alpha2 * h / beta2;
mat_A(l,l-1) = mat_A(l,l-2) + ci_x(x(l));

y = [zeros(l,1)];
mat_A
vec_b
% Solves tridigonal equations using Thomas Algorithm
y = thomas_algorithm(mat_A, vec_b);


plot(x,y)
legend('Solution for y = f(x)')
xlabel('X');
ylabel('Y');