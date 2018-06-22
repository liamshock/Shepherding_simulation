function createfile(x, y, v, u, filename)
    %creates a 1D of 2D data file with all the positions and velocities
    %for each time step.
    
    fid = fopen(filename, 'w'); %open file for writing
    
    timesteps = length(x(1,:));
    N = length(x(:,1));
    
    for t=1:timesteps-1 %print each timestep
        fprintf(fid, '\nTimestap %4.0f of %4.0f',t-1,timesteps-1);
        if y == zeros(size(y)) %print a 1D file
            fprintf(fid, '\n   N       X       V\n');
            for i=1:N
                fprintf(fid, '%4.0f %7.3f %7.3f\n', i, x(i,t), v(i,t));
            end
        else %print a 2D file
            fprintf(fid, '\n   N       X       Y       V       U\n');
            for i=1:N
                fprintf(fid, '%4.0f %7.3f %7.3f %7.3f %7.3f \n', i, x(i,t), y(i,t), v(i,t), u(i,t));
            end
        end
    end
    
    fclose(fid);

end