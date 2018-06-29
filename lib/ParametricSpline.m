function [ xt, yt, tt ] = ParametricSpline(x,y)
    n = length(x);
    t = zeros(n, 1);
    for i=2:n
        arc_length = sqrt((x(i)-x(i-1))^2 + (y(i)-y(i-1))^2);
        t(i) = t(i-1) + arc_length;
    end
    tt = linspace(0,1,1000);
    t=t./t(length(t));
    xt = spline(t, x, tt);
    yt = spline(t, y, tt);
end 
