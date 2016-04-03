clear; clc; close all 

% Draw a sphere
sphere
% Make the current axis box square in size
axis('square') 

% Define title and labels for reference
title('Rotation of a sphere...')
xlabel('x'); ylabel('y'); zlabel('z') 

% Modify azimuth (horizontal rotation) and update drawing
for az = -50 : .2 : 30
    view(az, 40)
    drawnow
end 

% Modify elevation (vertical rotation) and update drawing
for el = 40 : -.2 : -30
    view(30, el)
    drawnow
end 