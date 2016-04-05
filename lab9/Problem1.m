% Solve Possions Equation using Line by Line Algorithm
% 
% Author: Mudit Bachhawat
% Roll: 13MA20023
% Creation Date: 4th April, 2016
% Last Updated on: 4th April, 2016
% 
% 
% d2u/dx2 + d2u/dy2 = x^2 + y^2
% u(x ,y) = 0 for 0 <= x,y <= 1
% h = 1/4

d = 0.25;
dx = d;
dy = d;

x0 = 0;
y0 = 0;

xn = 1;
yn = 1;

x = x0:dx:xn;
y = y0:dy:yn;

i = length(x) - 1;
j = length(y) - 1;

% u(x,y) = 0 at the boundary
u = zeros(i+1, j+1, 1);

% Solution

ai_x = @(x)( 16 );
bi_x = @(x)( -64 );
ci_x = @(x)( 16 );
di_x = @(x,y,um1,u,up1)( x*x + y*y - 16*um1 - 16*up1 );

temp = u+1;
while max(max(abs(temp-u)))>1e-5,
    temp = u;
    for a = 2:i,

        mat_A = zeros(j-1, j-1);
        vec_b = zeros(j-1,1);

        % Preparing trigonal matrix
        for g=1:(j),
            mat_A(g,g) = bi_x(NaN);
            mat_A(g,g+1) = ci_x(NaN);
            mat_A(g+1,g) = ai_x(NaN);

            if (g~=j),
                vec_b(g) = di_x(x(a), y(g+1), u(a-1,g+1), u(a,g+1), u(a+1,g+1));
            end
        end


    %     Since u(x,y) = 0 on the boundary
    %     vec_b(1) = vec_b(1) - r*u0;
    %     vec_b(j-1) = vec_b(j-1) - r*un;

        mat_A = mat_A(1:j-1, 1:j-1);

        vec_b = vec_b(1:j-1);

        u(a, 2:j) = thomas_algorithm(mat_A, vec_b);


    end
end

mesh(x,y,u')
xlabel('X')
ylabel('Y')
zlabel('U(x)')