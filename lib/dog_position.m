function [x_dog_new, y_dog_new] = dog_position(x_bar, y_bar, phi, N_dogs, dog_dist)
% find the position of the dog that would move the herd in the correct
% direction
dx = zeros(N_dogs,1);
dy = zeros(N_dogs,1);
x_dog_new = zeros(N_dogs,1);
y_dog_new = zeros(N_dogs,1);

for k = 1:N_dogs
    dx(k) = dog_dist(k)*cos(phi);
    dy(k) = dog_dist(k)*sin(phi);
    x_dog_new(k)= x_bar - dx(k);
    y_dog_new(k) = y_bar - dy(k);
end

end

