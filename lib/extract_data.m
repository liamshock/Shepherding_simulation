function [A, B, Orient, Eccen, X_bar, Y_bar, Area] = extract_data(x, y, time, dt, N)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here


% Introductory message
fprintf('\n\n')
fprintf('--- Extracting the data ---')  
fprintf('\n')

% calculate the timesteps
timesteps = time/dt+1;

% preallocate arrays
A = zeros([timesteps 1]);
B = zeros([timesteps 1]);
Orient = zeros([timesteps 1]);
Eccen = zeros([timesteps 1]);
X_bar = zeros([timesteps 1]);
Y_bar = zeros([timesteps 1]);
Area = zeros([timesteps 1]);

% loop over timesteps
for t = 1:timesteps

    % take a snapshot of the positions at this time
    x_snapshot = x(:,t);
    y_snapshot = y(:,t);
    
    % get the convex hull area
    [~, hull_area] = convhull(x(:,t),y(:,t));
    Area(t) = hull_area;
    
    % Get the centre of mass variables
    Z = (1/N)*[sum(x(:,t)) sum(y(:,t))]'; % COM (\bar{x} \bar{y})
    X_bar(t) = Z(1);
    Y_bar(t) = Z(2);
    
    % fit ellipse
    ellipse_t = fit_ellipse(x_snapshot, y_snapshot);
    
    % get the orientation, eccentricity, major and minor axes
    A(t) = ellipse_t.long_axis;
    B(t) = ellipse_t.short_axis;
    Orient(t) = ellipse_t.phi;
    Eccen(t) = sqrt(1 - ((B(t)/A(t))^2));
        

end 

% end of the function
fprintf('Collective variables extracted')
fprintf('\n') 
end