function [error_modulus] = find_error_modulus(x_bar_init, y_bar_init, x_bar, y_bar, theta)
% find the modulus of the error vector


x_intersection = (y_bar_init - y_bar + x_bar*tan((pi/2) + theta) - x_bar_init*tan(theta)) / (tan((pi/2) + theta) - tan(theta));
y_intersection = x_intersection*tan(theta) + y_bar_init - x_bar_init*tan(theta);

error_modulus = norm([x_bar - x_intersection y_bar - y_intersection]);

end

