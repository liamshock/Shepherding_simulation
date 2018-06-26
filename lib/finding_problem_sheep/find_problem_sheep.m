function [problem_sheep_info, circle_centre] = find_problem_sheep(x_snapshot, y_snapshot)
% Get information about which sheep are deviating most for the circular
% shape
%
% -- Returns --
% (1) problem_sheep_info = a matrix (numSheepOnConvexHull, 3) 
%                          1st column = sheep x coord
%                          2nd column = sheep y coord
%                          3rd column = sheep distance to centre, with
%                                       descending distance among sheep
% (2) circle_centre = the centre of the 'inner circle' from the 
%                     'calculate_inner_circle' function
                

% So first we want to find the centre and radius of the 'inner circle'
[circle_centre, circle_radius] = calculate_inner_circle(x_snapshot, y_snapshot);

% get the points on the convex hull 
[K, ~] = convhull(x_snapshot, y_snapshot);

% get the number of points along the convex hull
numBoundaryPoints = size(K, 1);

% get a vector of the boundary points
boundaryPoints = [x_snapshot(K) y_snapshot(K)];

% find set of distances from the circle_centre to the boundary points
distances = zeros(numBoundaryPoints, 1);
for i = 1:numBoundaryPoints
    distances(i) = norm( boundaryPoints(i,:) - circle_centre );
end
    
% create an array to hold the output
problem_sheep_info = zeros(numBoundaryPoints, 3);
for i = 1:numBoundaryPoints
    problem_sheep_info(i,1) = boundaryPoints(i,1);
    problem_sheep_info(i,2) = boundaryPoints(i,2);
    problem_sheep_info(i,3) = distances(i);
end
    
% sort the list by distance i.e. "importance" - first row is most deviant
% sheep
problem_sheep_info = sortrows(problem_sheep_info, 3, 'descend');

end

