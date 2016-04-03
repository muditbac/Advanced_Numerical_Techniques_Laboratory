% Solve linear BVP and implement Spline Interpolation Method
% 
% 
% Author: Mudit Bachhawat
% Creation Date: 2nd April, 2016
% Last Updated on: 2nd April, 2016
% 
% 
% y'' + 2y' + y  = 30x
% y(0) = 0, y(1) = 0
% h = 1/2

h = 1/2;

x0 = 0;
xn = 1;
n = (xn - x0) / h;

x = x0:h:xn;

a = 2;
b = 1;
c = @(x)(30*x);

ax = @(x,wm1,w,wp1)([0 0; (-a/h) (a*h/6)]);
bx = @(x,wm1,w,wp1)([ b-a/h 1-h*a/3; b+a/h 1+h*a/3 ]);
cx = @(x,wm1,w,wp1)([a/h -h*a/6; 0 0]);

dx = @(x,wm1,w,wp1)([ c(x)  c(x) ]);

w = zeros(2, n+1);

w(1,1) = 0; % y0
w(1, end) = 0; % yn

mat_A = zeros(2,2,n-1, n-1);
vec_b = zeros(2,1, n-1);

hold on;


for i=2:n,
    index = i-1;
    mat_A(:,:,index,index) = bx(x(i), w(:,i-1), w(:,i), w(:,i+1));
    mat_A(:,:,index+1,index) = ax(x(i), w(:,i-1), w(:,i), w(:,i+1));
    mat_A(:,:,index,index+1) = cx(x(i), w(:,i-1), w(:,i), w(:,i+1));
    vec_b(:, index) = dx(x(i), w(:,i-1), w(:,i), w(:,i+1));
end

mat_A = mat_A(:,:,1:n-1,1:n-1);

p1 = 1-2*h/3;
mat_A(:,:,1,1) = mat_A(:,:,1,1) + (h*a/(6*p1))*[0 0; -2/h h/3];

p2 = 1+2*h/3;
mat_A(:,:,end,end) = mat_A(:,:,end,end) + (-h*a/(6*p2))*[2/h -h/3; 0 0];
vec_b(:,:,end) = vec_b(:,:,end) + (-h*a/(6*p2))*[0;c(x(end))];

w = thomas_algorithm_block(mat_A,vec_b);

% Using backward difference method
% w = w + dw;

M0 = (h*w(2,1)/3 - 2*w(1,1)/h)/p1;
Mn = (-h*w(2,end)/3 + 2*w(1,end)/h + c(x(end)))/p2;
% w(1,end) = yn;
w = [[0;M0] w [0;Mn]];


p = plot(x,w(1,:),'o-', 'LineWidth',1);


xx = x0:h/100:xn;
% yy = spline(x,w(1,:), xx);
yy = spline_interpolation(x, w(1,:), w(2,:), xx);

% set(p, 'LineWidth',1.5, 'LineStyle', '-')
plot(xx,yy);
xlabel('X Axis');
ylabel('Y Axis');
legend('Normal Values','Spline Interpolation')
