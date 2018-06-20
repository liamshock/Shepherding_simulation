%This script plots a bunch of particles and their behavior
%The script can read both 1D and 2D files

%%% DATA TO CHANGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = 'input2d.txt'; %name of the data file which should be loaded
sigma = 1;              %particle diameter
epsilon = 0.1;           %material parameter
L = 50;                 %attractive length in particle diameters
G = 5;                  %attractive strength
lambda = 0.1;           %damping coefficient
rep = -1;               %repulsion from shepherd (set -1 for no shepherd)
time = 50;              %total time
dt = 0.001;             %time step
method = 'Verlet';      %integration method
borders = [00 30 00 30];%set [x0 x1 y0 y1] to turn borders on (particles bounce), 0 is off
output = '';            %name of the output file. Set '' for no output
movie = 'output.mp4';   %movie output name. Set '' for no movie output
fps = 10;               %FPS for movie
ApplyBC = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


timesteps = time/dt+1; %number of timesteps to be calculated
[x y v u m variables N] = loadfile( filename, timesteps ); %load the data

V = zeros(N,timesteps); %allocate potential energy array
T = zeros(N,timesteps); %allocate kinetic energy array
f_x = zeros(1,N);       %allocate force matrix x-direction
f_y = zeros(1,N);       %allocate force matrix y-direction
V_j = zeros(1,N);       %allocate potential matrix

% create a position array for the sentinel

%%% Calcualtion part %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t=1:time/dt+1
    for i=1:N
        for j=1:N %set force/energy interactions for each particle
            if i~=j %don't take self inducting terms into account
                
                %compute the distance
                rij = sqrt((x(i,t)-x(j,t))^2+ (y(i,t)-y(j,t))^2);
                
                %compute unit normal vector, n(1)=n_x, n(2)=n_y
                n = [ x(j,t)-x(i,t) y(j,t)-y(i,t) ]/norm([ x(j,t)-x(i,t) y(j,t)-y(i,t) ]);
                
                %compute effective attractive strength
                Geff = G*(1+(j==1)*(-rep-1));
                
                %for x direction
               %f_x(j) = n(1)*-4*epsilon*((6*R*(L*sigma)^6)/rij^7 - (12*(L*sigma)^12)/rij^13);
                f_x(j) = n(1)*-4*epsilon*(exp(-rij/sigma)/sigma - Geff*exp(-rij/(L*sigma))/(L*sigma));
                
                %for y direction
               %f_y(j) = n(2)*-4*epsilon*((6*(L*sigma)^6)/rij^7 - (12*(L*sigma)^12)/rij^13);
                f_y(j) = n(2)*-4*epsilon*(exp(-rij/sigma)/sigma - Geff*exp(-rij/(L*sigma))/(L*sigma));
                
                %potential energy
               %V_j(j) = 4*epsilon*((sigma/rij)^12-(L*sigma/rij)^6);
                V_j(j) = 4*epsilon*(Geff*exp(-rij/(L*sigma)) - exp(-rij/sigma));
            
            else %if i=j: force=0, potential=0
                f_x(j) = 0;
                f_y(j) = 0;
                V_j(j) = 0;
            end
        end
        
        %Compute the position and velocity using Euler or Verlet
        if (strcmp('Euler',method) || t==1) %first timestap is always Euler   
            % Forward Euler x
            x(i,t+1) = x(i,t) + v(i,t)*dt;
            v(i,t+1) = v(i,t) + (sum(f_x)-lambda*v(i,t))/m(i)*dt;
            
            % Forward Euler y
            y(i,t+1) = y(i,t) + u(i,t)*dt;
            u(i,t+1) = u(i,t) + (sum(f_y)-lambda*v(i,t))/m(i)*dt;
            
        elseif(strcmp('Verlet',method))
            %Verlet algorithm x
            x(i,t+1) = -x(i,t-1) + 2*x(i,t) + (sum(f_x)-lambda*v(i,t))/m(i)*dt^2;
            v(i,t) = (x(i,t-1) - x(i,t+1))/(2*dt);
           
            %Verlet algorithm y
            y(i,t+1) = -y(i,t-1) + 2*y(i,t) + (sum(f_y)-lambda*v(i,t))/m(i)*dt^2;
            u(i,t) = (y(i,t-1) - y(i,t+1))/(2*dt);
        else
            disp('No valid integration method given!');
        end
        
        %%% If borders are on %%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        if ApplyBC
            if (x(i,t+1) < borders(1)+sigma/2) || (x(i,t+1) > borders(2)-sigma/2) %x
                x(i,t+1) = 2*x(i,t) - x(i,t+1);
                if (strcmp('Euler',method))
                    v(i,t+1) = -v(i,t+1);
                elseif(strcmp('Verlet',method))
                    v(i,t) = -v(i,t);
                end     
            end
            if (y(i,t+1) < borders(3)+sigma/2) || (y(i,t+1) > borders(4)-sigma/2) %y
                y(i,t+1) = 2*y(i,t) - y(i,t+1);
                if (strcmp('Euler',method))
                    u(i,t+1) = -u(i,t+1);
                elseif(strcmp('Verlet',method))
                    u(i,t) = -u(i,t);
                end   
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %Energy calculations
        V(i,t) = sum(V_j); %Potential energy: sum of influence of other particles
        T(i,t) = 1/2*m(i)*(v(i,t)^2+u(i,t)^2); %Kinetic energy: 1/2mv^2
    end
end

%%% Calculate and plot energy %%%%%%%%%%%%%%%%%%%%%%%%
V_sum = sum(V)/2; %sum the potential energy of all particles times a half since matrix is mirrored
T_sum = sum(T) ;%1/2.*m'.*(sum(u.^2)+sum(v.^2)); %sum from all particles: 0.5*m*v^2

fig=figure(1);
hold on
plot(0:dt:time,T_sum+V_sum)
plot(0:dt:time,T_sum,'r')
plot(0:dt:time,V_sum,'k')
title('The total energy of the system','FontSize',14);
xlabel('Time','FontSize',12);
ylabel('Energy','FontSize',12);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Create file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~strcmp('', output) %check if output file is provided
    createfile(x, y, v, u, output);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Set the axis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ApplyBC %if borders are given set the axis equal to the borders
    ax_x = borders(1:2);
    ax_y = borders(3:4);
else %if not set the axis automatically to the highest traveled distance of the particles
    ax_x = [min(min(x))-0.5*sigma max(max(x))+0.5*sigma]; %x axis
    if variables == 3
        ax_y = [-1 1]; %y axis for 1D
    else
        ax_y = [min(min(y))-0.5*sigma max(max(y))+0.5*sigma]; %y axis for 2D
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Draw the circles and make the movie%%%%%%%%%%%%%%%
if ~(strcmp('',movie))
    fig=figure(2);
    mov = VideoWriter(movie, 'MPEG-4');

    for k=1:1/dt/fps:timesteps %pick the points for the correct number of fps
        circle = linspace(0,2*pi,100); %create a circle
        hold on;
        for n=1:N %plot N circles...
            xx = sigma/2*cos(circle)+x(n,round(k)); %x-location of particle                  
            yy = sigma/2*sin(circle)+y(n,round(k)); %y-location of particle

            if variables == 3
                color = [abs(v(n,round(k)))/max(max(abs(v))) 0 0];
            else
                color = [abs(v(n,round(k)))/max(max(abs(v))) 0 abs(u(n,round(k)))/max(max(abs(u)))];
            end
            %color is based on the max velocity in v and u direction
            %red = high v, blue = high u, or combo

            fill(xx,yy,color);
        end   
        
        % get snapshots
        x_snapshot = x(:,k);
        y_snapshot = y(:,k);
        
        % get the convex hull points
        K = convhull(x_snapshot, y_snapshot);
        
        % plot the convex hull
        plot(x_snapshot(K), y_snapshot(K),'b','linewidth',3);
        %hold on;
        
        % fit ellipse
        if k ~= 1
            position_snapshot = squeeze(cat(3,x_snapshot(K), y_snapshot(K)))';
            [z, a, b, alpha] = fitellipse(position_snapshot);
            plotellipse(z, a, b, alpha)
        end

        set(gca,'DataAspectRatio',[1 1 1]); %set aspect ratio x:y to 1:1
        axis([ax_x ax_y]); %set the axis
        %cleahold on;

        %make the movie
        F=getframe(fig);
        open(mov);
        writeVideo(mov,F);

        %clear the figure for removing traces
        clf(fig);
    end

    %close the movie when done
    %mov = close(mov);
end
close(mov);


%%%%%%% Create the position and velocity arrays %%%%%%%
position = cat(3,x,y);
velocity = cat(3,v,u);

%%%%%%% Plot the boundary of the ellipse %%%%%%

x_snapshot = squeeze(x(:,2000));
y_snapshot = squeeze(y(:,2000));

%fig = figure(3);
%plot(x_snapshot, y_snapshot, '.');
%k = boundary(x_snapshot,y_snapshot);
%hold on;
%plot(x_snapshot(k),y_snapshot(k));


%%%%%% Fit a minimum volume ellipse to the data %%%%%
% position_snapshot = squeeze(position(:,2000,:))';
% [A, c] = MinVolEllipse(position_snapshot, 0.01);
% 
% plot(x_snapshot, y_snapshot, '*');
% hold on;
% Ellipse_plot(A,c);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
