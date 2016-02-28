% Solve non linear BVP using Newton Linearization Technique
% 
% 
% Author: Mudit Bachhawat
% Creation Date: 25th Feburary, 2015
% Last Updated on: 25th Feburary, 2015
% 
% 
% y'' = .5*(1+x+y)^3
% y(0) = 0, y(1) = 0
% 

h = .02;

a = @(x,ym1,y,yp1)(1);
b = @(x,ym1,y,yp1)(-2-1.5*h*h*(1+x+y)^2);
c = @(x,ym1,y,yp1)(1);

d = @(x,ym1,y,yp1)(h/2*h*(1+x+y)^3-yp1+2*y-ym1);

x = 0:h:1;
n = length(x)-1;

f = @(x)(x*(1-x));

mat_A = zeros(n-1, n-1);
vec_b = zeros(1, n-1);
y = zeros(1, n+1);

for i=1:n+1,
    y(i) = f(x(i));
end

legend_vals = {'k =0'};
plot(x,y,'-.', 'LineWidth',1);
hold on;
% y = y(1:n+1);

% Random initialization
dy = 1;
k = 0;
while max(abs(dy))>1e-5,
    k = k+1;

    for i=2:n,
        index = i-1;
        mat_A(index,index) = b(x(i), y(i-1), y(i), y(i+1));
        mat_A(index+1,index) = a(x(i), y(i-1), y(i), y(i+1));
        mat_A(index,index+1) = c(x(i), y(i-1), y(i), y(i+1));
        vec_b(index) = d(x(i), y(i-1), y(i), y(i+1));
    end

    mat_A = mat_A(1:n-1,1:n-1);

%     dy = vec_b/mat_A;
    dy = thomas_algorithm(mat_A,vec_b');
    dy = [0 dy' 0];
    y = y + dy;

    
    p = plot(x,y,'-.', 'LineWidth',1);
    
    legend_vals = [legend_vals (strcat('k = ',num2str(k)))];

end

set(p, 'LineWidth',1.5, 'LineStyle', '-')
xlabel('X Axis');
ylabel('Y Axis');


legend(legend_vals)

figure
plot(x,y)
legend('y = f(x)')