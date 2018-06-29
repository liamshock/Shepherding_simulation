function [] = make_movie(outpath, dt, time, fps, N, ApplyBC, borders, x, y, u, v, x_dog, ...
                         y_dog, u_dog, v_dog, X_T, Y_T, N_dogs, ...
                         DBeta, tol,Wx,Wy, d1, d2, d3, d4)
% Make and write a movie of the simulation results


% Introductory message
fprintf('\n\n')
fprintf('--- Making the movie ---')  
fprintf('\n')

% set the output name of the movie
movieName = 'movie.mp4';    
moviePath = fullfile(outpath, movieName);

% recalculate the timesteps
timesteps = time/dt+1; %number of timesteps to be calculated

% set the axes 
% if ApplyBC %if borders are given set the axis equal to the borders
%     ax_x = borders(1:2);
%     ax_y = borders(3:4);
% else %if not set the axis automatically to the highest traveled distance of the particles
%     ax_x = [min([min(min(x)),borders(1),min(x_dog)])-0.5 max([max(max(x)),borders(2),max(x_dog)])+0.5]; %x axis
%     ax_y = [min([min(min(y)),borders(3),min(y_dog)])-0.5 max([max(max(y)),borders(4),max(y_dog)])+0.5]; %y axis for 2D
% end



% Draw the circles and make the movie
%  fig=figure('visible', 'off');

mov = VideoWriter(moviePath, 'MPEG-4');
mov.FrameRate = fps;

% get the total number of frames we want to write
numFrames = fps*time;

% get the interval that we want for the writing frames
interval = timesteps/numFrames;

% set initial target index
target_idx = 1;
x_T = X_T(target_idx);
y_T = Y_T(target_idx);

% write the movie, keep tracking of frames written
frames_written = 0;
for k=1:interval:timesteps 
figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off')    
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
    %plot N circles
    for n=1:N
        THETHA = atan2(u(n,round(k)),v(n,round(k)));
        xx = 0.12*cos(circle)*cos(THETHA) - 0.12*sin(circle)*sin(THETHA) + x(n,round(k)); %x-location of particle                  
        yy = 0.12*cos(circle)*sin(THETHA) + 0.12*sin(circle)*cos(THETHA) + y(n,round(k)); %y-location of particle
        %color = [abs(v(n,round(k)))/max(max(abs(v))) 0 abs(u(n,round(k)))/max(max(abs(u)))];
        %color is based on the max velocity in v and u direction
        %red = high v, blue = high u, or combo

        fill(xx,yy,'w');
        hold on 
    end 

    hold on
    % get snapshots
    x_snapshot = x(:,round(k));
    y_snapshot = y(:,round(k));

%     % get the convex hull points
%     K = convhull(x_snapshot, y_snapshot);
% 
%     % plot the convex hull
%     plot(x_snapshot(K), y_snapshot(K),'b','linewidth',3);
%     %hold on;

    % plot the predator
    for i = 1:N_dogs
        THETHA_z = atan2(u_dog(i,round(k)), v_dog(i,round(k)));
        xx_z = 0.15*cos(circle)*cos(THETHA_z) - 0.15*sin(circle)*sin(THETHA_z) + x_dog(i,round(k)); %x-location of particle                  
        yy_z = 0.15*cos(circle)*sin(THETHA_z) + 0.15*sin(circle)*cos(THETHA_z) + y_dog(i,round(k)); %y-location of particle
        fill(xx_z,yy_z,'k');
    end
        %hold on

    % fit ellipse and plot 
    % NOTE: We have already fitted the ellipse in the 'extract_data'
    % function. We re-do it here because it has the ability to add itself
    % to a plot. Inefficient, but it works.
%     ellipse_t = fit_ellipse(x_snapshot, y_snapshot, fig);

    % set the axes
    %set(gca,'DataAspectRatio',[1 1 1]); %set aspect ratio x:y to 1:1
    %axis([ax_x ax_y]); %set the axis
    %cleahold on;
    %title('\fontsize{14} Circular motion of dog');
%     legend([p1,p2,p3],'Sheep','Convex Hull','Dog', 'Location','northwest')
    
    % plot the trajectory line
    %plot([x_bar_init x_T], [y_bar_init y_T]);
    
    
    
    
    
 
    
    % plot check
    %x_curr = x(:,1); y_curr = y(:,1); xd_curr = x_dog(:,1); yd_curr = y_dog(:,1);
    %     figure 
    %     plot(x_curr,y_curr,'*k'); axis equal tight;
    %     hold on
    %     plot(xd_curr,yd_curr,'or');
    %     K  = convhull(x_curr,y_curr);
    %     plot(x_curr(K),y_curr(K),'b'); 
    xbar = mean(x_snapshot); ybar = mean(y_snapshot); 
    
                 % decide on targets
         if sqrt((x_T - xbar)^2 + (y_T - ybar)^2) < tol
             % dont go beyond number of targets
             target_idx = min(target_idx + 1, size(X_T, 1));
             x_T = X_T(target_idx);
             y_T = Y_T(target_idx);
         end 
         
    %     plot(xbar,ybar,'og'); 
    %     plot([xbar,x_T],[ybar,y_T],'g'); 
    Beta = atan2((-ybar+y_T),(-xbar+x_T));
    Edel = (2*pi-(2*DBeta))/4;
    Th1 = Beta+DBeta+Edel;
    Th2 = Th1+Edel;
    Th3 = Th2+Edel;
    maxDis = max(sqrt((x_snapshot-xbar).^2+(y_snapshot-ybar).^2));
    radDogDist = 1.1*maxDis;
    px1 = xbar + radDogDist*cos(Th1); py1 = ybar + radDogDist*sin(Th1);
    px2 = xbar + radDogDist*cos(Th2); py2 = ybar + radDogDist*sin(Th2);
    px3 = xbar + radDogDist*cos(Th3); py3 = ybar + radDogDist*sin(Th3);
    %targetXdog = [px1,px2,px3]; targetYdog = [py1,py2,py3];
%     plot([xbar,px1],[ybar,py1],'x-r'); 
%     plot([xbar,px2],[ybar,py2],'x-r'); 
%     plot([xbar,px3],[ybar,py3],'x-r'); 
        
        
        
        
        
        
        
        % targeted motion 
%     Wx = 50; Wy = 30; d1 = 10; d2 = 12; d3 = 2; d4 = 6; 
    %Fullborders = [0 Wx 0 Wy];                      
    gpx1 = Wx-d2; gpx2 = gpx1+0.5*(d2-d4); gpx3 = gpx2+d4; gpx4 = Wx; gpx5 = gpx1; gpx6 = gpx1-d3;
    gpy1 = Wy-d1; gpy2 = gpy1; gpy3 = gpy1; gpy4 = gpy1; gpy5 = Wy; gpy6 = gpy5;
    %x_T = 0.5*(gpx1+gpx2);                                 % desired x-comp of flock COM
    %y_T = 0.5*(gpy1+gpy6);                                % desired x-comp of flock COM

    % inspection figure 
    %figure 
    plot([gpx1,gpx6],[gpy1,gpy6],'k','linewidth',3);%hold on
    plot([gpx1,gpx2],[gpy1,gpy2],'k','linewidth',3);
    plot([gpx3,gpx4],[gpy3,gpy4],'k','linewidth',3);
    plot([0,Wx],[0,0],'k','linewidth',3);
    plot([Wx,Wx],[0,Wy],'k','linewidth',3);
    plot([Wx,0],[Wy,Wy],'k','linewidth',3);
    plot([0,0],[Wy,0],'k','linewidth',3);
    axis equal tight 
    %hold on
%     plot([borders(1:2)],[borders(3),borders(3)],'r','linewidth',2);
%     plot([borders(2),borders(2)],[borders(3:4)],'r','linewidth',2);
%     plot([borders(2),borders(1)],[borders(4),borders(4)],'r','linewidth',2);
%     plot([borders(1),borders(1)],[borders(4),borders(3)],'r','linewidth',2);

    xpthIc = mean(borders(1:2));
    ypthIc = mean(borders(3:4));
    plot(xpthIc,ypthIc,'*r');
    plot(X_T(end),Y_T(end),'*r');
    
    % plot the trajectory
    [xt, yt, tt ] = ParametricSpline(X_T, Y_T);
    plot(xt, yt, '--r','linewidth',2);
    %plot(X_T(1:end-1),Y_T(1:end-1),'--r','linewidth',2);
    
    
    % display the current target
    plot(X_T, Y_T, 'go');
    set(gca,'fontsize',18);
    set(gcf,'color','w');
    set(gca,'color', [0.6, 0.9, 0.55]);



    %make the movie
    F=getframe(gcf);
    open(mov);
    writeVideo(mov,F);

    %clear the figure for removing traces
    close all
end
%close the movie when done
close(mov)



% end of the function
fprintf('Finished writing the movie')
fprintf('\n')
end

