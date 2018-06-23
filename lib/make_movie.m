function [] = make_movie(outpath, sigma, dt, time, fps, N, ApplyBC, borders, x, y, u, v, x_dog, y_dog, u_dog, v_dog)
% Make and write a movie of the simulation results


% Introductory message
fprintf('\n\n')
fprintf('--- Making the movie---')  
fprintf('\n')

% set the output name of the movie
movieName = 'movie.mp4';    
moviePath = fullfile(outpath, movieName);

% recalculate the timesteps
timesteps = time/dt+1; %number of timesteps to be calculated

% set the axes 
if ApplyBC %if borders are given set the axis equal to the borders
    ax_x = borders(1:2);
    ax_y = borders(3:4);
else %if not set the axis automatically to the highest traveled distance of the particles
    ax_x = [min(min(x))-0.5*sigma max(max(x))+0.5*sigma]; %x axis
    ax_y = [min(min(y))-0.5*sigma max(max(y))+0.5*sigma]; %y axis for 2D
end



% Draw the circles and make the movie
fig=figure('visible', 'off');
mov = VideoWriter(moviePath, 'MPEG-4');
mov.FrameRate = fps;

% get the total number of frames we want to write
numFrames = fps*time;

% get the interval that we want for the writing frames
interval = timesteps/numFrames;

% write the movie, keep tracking of frames written
frames_written = 0;
for k=1:interval:timesteps 
    
    % update the number of frames written
    frames_written = frames_written + 1;
    
    % print quarterly updates
    if frames_written == round(numFrames/10)
        fprintf('1/10 complete \n')
    elseif frames_written == round(numFrames/5)
        fprintf('2/10 complete \n')
    elseif frames_written == round( (3/10)*numFrames)
        fprintf('3/10 complete \n')
    elseif frames_written == round( (4/10)*numFrames)
        fprintf('4/10 complete \n')   
    elseif frames_written == round(numFrames/2)
        fprintf('5/10 complete \n')   
    elseif frames_written == round( (6/10)*numFrames)
        fprintf('6/10 complete \n')  
    elseif frames_written == round( (7/10)*numFrames)
        fprintf('7/10 complete \n')
    elseif frames_written == round( (8/10)*numFrames)
        fprintf('8/10 complete \n')
    elseif frames_written == round( (9/10)*numFrames)
        fprintf('9/10 complete \n')
    end
    
    circle = linspace(0,2*pi,100); %create a circle
    hold on;
    %plot N circles
    for n=1:N
        THETHA = atan2(u(n,round(k)),v(n,round(k)));
        xx = sigma/2*cos(circle)*cos(THETHA) - sigma/4*sin(circle)*sin(THETHA) + x(n,round(k)); %x-location of particle                  
        yy = sigma/2*cos(circle)*sin(THETHA) + sigma/4*sin(circle)*cos(THETHA) + y(n,round(k)); %y-location of particle
        color = [abs(v(n,round(k)))/max(max(abs(v))) 0 abs(u(n,round(k)))/max(max(abs(u)))];
        %color is based on the max velocity in v and u direction
        %red = high v, blue = high u, or combo

        p1 = fill(xx,yy,color);
    end 

    hold on
    % get snapshots
    x_snapshot = x(:,round(k));
    y_snapshot = y(:,round(k));

    % get the convex hull points
    K = convhull(x_snapshot, y_snapshot);

    % plot the convex hull
    p2 = plot(x_snapshot(K), y_snapshot(K),'b','linewidth',3);
    hold on;

    % plot the predator
    THETHA_z = atan2(u_dog(round(k)), v_dog(round(k)));
    xx_z = sigma/2*cos(circle)*cos(THETHA_z) - sigma/4*sin(circle)*sin(THETHA_z) + x_dog(round(k)); %x-location of particle                  
    yy_z = sigma/2*cos(circle)*sin(THETHA_z) + sigma/4*sin(circle)*cos(THETHA_z) + y_dog(round(k)); %y-location of particle
    p3 = fill(xx_z,yy_z,'c');
    hold on

    % fit ellipse and plot 
    % NOTE: We have already fitted the ellipse in the 'extract_data'
    % function. We re-do it here because it has the ability to add itself
    % to a plot. Inefficient, but it works.
    ellipse_t = fit_ellipse(x_snapshot, y_snapshot, fig);

    % set the axes
    set(gca,'DataAspectRatio',[1 1 1]); %set aspect ratio x:y to 1:1
    axis([ax_x ax_y]); %set the axis
    %cleahold on;
    %title('\fontsize{14} Circular motion of dog');
    legend([p1,p2,p3],'Sheeps','convex hull','Dog')

    %make the movie
    F=getframe(fig);
    open(mov);
    writeVideo(mov,F);

    %clear the figure for removing traces
    clf(fig);
end
%close the movie when done
close(mov)



% end of the function
fprintf('Finished writing the movie')
fprintf('\n')
end

