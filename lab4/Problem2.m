% Solve non linear BVP using Newton Linearization Technique
% 
% 
% Author: Mudit Bachhawat
% Creation Date: 25th Feburary, 2015
% Last Updated on: 27th Feburary, 2015
% 
% 
% y''' + yy'' + 1 - (y')^2 = 0
% y(0) = 0, y'(0) = 0
% y'(10) = 1

h = 0.1;

a = @(x,wm1,w,wp1)([0 (1/(h*h)-w(1)/(2*h)); (-1) (-h/2)]);
b = @(x,wm1,w,wp1)([ (wp1(2)-wm1(2))/(2*h) (-2/(h*h)-2*w(2)); (1) (-h/2) ]);
c = @(x,wm1,w,wp1)([0 (1/(h*h)+w(1)/(2*h)); 0 0]);

d = @(x,wm1,w,wp1)([(-1+(w(2))^2-w(1)*(wp1(2)-wm1(2))/(2*h)-(wp1(2)+wm1(2)-2*w(2))/(h*h))
    (-w(1)+wm1(1)+(h/2)*(w(2)+wm1(2)))]);


x = 0:h:10;
n = length(x)-1;

f = @(x)(5*x^2 - x^3/3 + x^2/20);
F = @(x)(x*(10-x)+x/10);

w = zeros(2, n+1);

for i=1:n+1,
    w(:,i) = [f(x(i)) F(x(i))];
end

mat_A = zeros(2,2,n-1, n-1);
vec_b = zeros(2,1, n-1);

legend_vals = {'k =0'};
plot(x,w(1,:),'-.', 'LineWidth',1);
hold on;

dw = 1;
k = 0;
while max(max(abs(dw)))>1e-5,
    k = k+1;

    for i=2:n,
        index = i-1;
        mat_A(:,:,index,index) = b(x(i), w(:,i-1), w(:,i), w(:,i+1));
        mat_A(:,:,index+1,index) = a(x(i), w(:,i-1), w(:,i), w(:,i+1));
        mat_A(:,:,index,index+1) = c(x(i), w(:,i-1), w(:,i), w(:,i+1));
        vec_b(:, index) = d(x(i), w(:,i-1), w(:,i), w(:,i+1));
    end

    mat_A = mat_A(:,:,1:n-1,1:n-1);

    dw = thomas_algorithm_block(mat_A,vec_b);
    
    % Using backward difference method
    dw = [[0;0] dw [0;0]];
    w = w + dw;

    yn = (2*h*w(2,end)+4*w(1,end-1)-w(1,end-2))/3;
    w(1,end) = yn;

    
    p = plot(x,w(1,:),'-.', 'LineWidth',1);
    
    legend_vals = [legend_vals (strcat('k = ',num2str(k)))];    
end

set(p, 'LineWidth',1.5, 'LineStyle', '-')
xlabel('X Axis');
ylabel('Y Axis');


legend(legend_vals)

figure
plot(x,w)
legend('y = f(x)','y = f''(x)')