function [x_dog_new, y_dog_new] = dog_position(x_bar, y_bar, phi, dog_dist)
% find the position of the dog that would move the herd in the correct
% direction

dx = dog_dist*cos(phi);
dy = dog_dist*sin(phi);
x_dog_new = x_bar - dx;
y_dog_new = y_bar - dy;
end

