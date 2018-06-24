function [x, y, v, u, m, variables, N] = loadfile( filename , timesteps)
    %load a file and returns the the position and velocity data
    %also prelocates the memory by setting all following timesteps to 0

    file = load( filename );
    variables = size(file,2); %number of variables
    N = size(file,1);    %number of particles

    %prelocate memory by setting all variables to zero
    x = zeros(N, timesteps);
    y = zeros(N, timesteps);
    v = zeros(N, timesteps);
    u = zeros(N, timesteps);
    
    %determine 1D or 2D situation and set the start conditions
    if variables == 3        %1D
        x(:,1) = file(:,1);
        v(:,1) = file(:,2);
        m      = file(:,3);
    elseif variables == 5    %2D
        x(:,1) = file(:,1);
        y(:,1) = file(:,2);
        v(:,1) = file(:,3);
        u(:,1) = file(:,4);
        m      = file(:,5);
    end
end