function [a, b] = extract_ellipse_parameters(A, C)
% extract_ellipse_parameters: return the ellipse information
% 
% This function should be used to return the ellipse information from the
% ellipse fitting by MinVolEllipse


% "singular value decomposition" to extract the orientation and the
% axes of the ellipsoid
[U D V] = svd(A);

% get the major and minor axes
a = 1/sqrt(D(1,1));
b = 1/sqrt(D(2,2));

end

