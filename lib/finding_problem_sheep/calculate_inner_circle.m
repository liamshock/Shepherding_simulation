function [circle_centre, circle_radius] = calculate_inner_circle(x_snapshot, y_snapshot)
% find the largest circle that can fit entirely within the convex hull
%
% ---- Returns ----
% circle_centre -> a (1,2) vector holding the centre of the circle
% circle_radius -> a scalar, the radius of the circle


% get the points on the convex hull 
[K, ~] = convhull(x_snapshot, y_snapshot);

% get the number of points along the convex hull
numBoundaryPoints = size(K, 1);

% get a vector of the boundary points
boundaryPoints = [x_snapshot(K) y_snapshot(K)];

% get centre of the circle, which is the COM of the points on the convex hull
circle_centre = [ (1/numBoundaryPoints)*sum(x_snapshot(K)) (1/numBoundaryPoints)*sum(y_snapshot(K)) ];

% find set of distances from the COM to the boundary points
distances = zeros(numBoundaryPoints, 1);
for i = 1:numBoundaryPoints
    distances(i) = norm( boundaryPoints(i,:) - circle_centre );
    
% the radius of the circle will be the smallest distance
circle_radius = min(distances);
end

