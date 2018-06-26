function plot_inner_circle_and_convex_hull(x_snapshot, y_snapshot)
% This function is purely for illustration
% It is used to test the 'calculate_inner_circle' function


% ----------calculations ---------

% calculate the inner circle
[circle_centre, circle_radius] = calculate_inner_circle(x_snapshot, y_snapshot);

% calculate the convex hull
K = convhull(x_snapshot, y_snapshot);



% ---------plotting ----------------

% set a figure
fig = figure(1);

% plot the points
p1 = plot(x_snapshot, y_snapshot, '*b');
hold on

% plot the convex hull
p2 = plot(x_snapshot(K), y_snapshot(K),'b','linewidth',3);
%hold on;

% plot the cricle
ang=0:0.01:2*pi; 
xp=circle_radius*cos(ang);
yp=circle_radius*sin(ang);
plot(circle_centre(1)+xp, circle_centre(2)+yp);
hold off
end