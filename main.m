%This script plots a bunch of particles and their behavior
%The script can read both 1D and 2D files

close all
clc

%%% DATA TO CHANGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set a path to save all of the results from the simulation
outpath = '/Users/liam/output_parfor';
if ~exist(outpath, 'dir')
  mkdir(outpath);
end

% set the output name of the movie
movieName = 'movie.mp4';    
moviePath = fullfile(outpath, movieName);

rng(4);
sigma = 1;                 % particle diameter
epsilon = 15;              % material parameter
L = 20;                    % attractive length in particle diameters
l = 5;                     % repulsive length from the predator
G = 0.1;                   % attractive strength
lambda = 1;                % damping coefficient
c = 0.4;                   % predator-prey repulsion
%rep = -1;                 % repulsion from shepherd (set -1 for no shepherd)
time = 1;                 % total time
dt = 0.001;                % time step
method = 'Kinematic';      % integration method
borders = [-30 30 -30 30]; % set [x0 x1 y0 y1] to turn borders on (particles bounce), 0 is off
fps = 10;                  % FPS for movie
ApplyBC = false;

N = 20; % num of particles 
m = 0.1*ones(N,1); % mass of the particles

timesteps = time/dt+1; %number of timesteps to be calculated

% initialize variables 
x = zeros(N,timesteps); y = x;  u = x; v = x;

% Init Position and velocity of a single Predator
thet0_z = linspace(0,4*pi,timesteps); %2*pi*rand(1);
x_z = 5.*cos(thet0_z)+15; y_z = 5.*sin(thet0_z)+15;
v_z = -5.*sin(thet0_z); u_z = 5.*cos(thet0_z);

% Initial conditions and velocities 
x0 = borders(1) + (borders(2)-borders(1))*rand(N,1); % x - initial positions of agents form uniform random distr 
y0 = borders(3) + (borders(4)-borders(3))*rand(N,1); % y - initial positions of agents form uniform random distr 

% initial velocities 
thet0 = 2*pi*rand(N,1); 
spd0 = rand(N,1); 
u0 = spd0.*cos(thet0);
v0 = spd0.*sin(thet0);
variables = 5;

% place the Ic in their corresponding matrices 
x(:,1) = x0;
y(:,1) = y0;
u(:,1) = u0;
v(:,1) = v0;

% initialize more variables
V = zeros(N,timesteps); %allocate potential energy array
T = zeros(N,timesteps); %allocate kinetic energy array
f_x = zeros(1,N);       %allocate force matrix x-direction
f_y = zeros(1,N);       %allocate force matrix y-direction
V_j = zeros(1,N);       %allocate potential matrix

A = zeros([timesteps 1]);
B = zeros([timesteps 1]);
save_Z = zeros(timesteps,2);
save_X_bar = zeros([timesteps 1]);
save_Y_bar = zeros([timesteps 1]);
Alpha = zeros([timesteps 1]);
Orient = zeros([timesteps 1]);
Eccen = zeros([timesteps 1]);
Z  = zeros(2,timesteps);
X_bar = zeros([timesteps 1]);
Y_bar = zeros([timesteps 1]);
Area = zeros([timesteps 1]);




%%%%%%%%%%%%%%%%%%% Calcualtion part %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for t=1:time/dt+1
    for i=1:N
        for j=1:N %set force/energy interactions for each particle
            if i~=j %don't take self inducting terms into account
                
                %compute the distance
                rij = sqrt((x(i,t)-x(j,t))^2+ (y(i,t)-y(j,t))^2);
                riz = sqrt((x(i,t)-x_z(t))^2+(y(i,t)-y_z(t))^2);
                
                %compute unit normal vector, n(1)=n_x, n(2)=n_y
                n = [ x(j,t)-x(i,t) y(j,t)-y(i,t) ]/norm([ x(j,t)-x(i,t) y(j,t)-y(i,t) ]);
                n_z = [ x_z(t)-x(i,t) y_z(t)-y(i,t) ]/norm([ x_z(t)-x(i,t) y_z(t)-y(i,t) ]);
                
                %for x direction
                f_x(j) = -4*epsilon*(n(1)*exp(-rij/sigma)/sigma - n(1)*G*exp(-rij/(L*sigma))/(L*sigma) + n_z(1)*c*exp(-riz/(l*sigma))/(l*sigma));
                
                %for y direction
                f_y(j) = -4*epsilon*(n(2)*exp(-rij/sigma)/sigma - n(2)*G*exp(-rij/(L*sigma))/(L*sigma) + n_z(2)*c*exp(-riz/(l*sigma))/(l*sigma));
                
                %potential energy
                V_j(j) = 4*epsilon*(-sigma*exp(-rij/sigma) + G*exp(-rij/(L*sigma)) - c*exp(-riz/(l*sigma)));
            
            else %if i=j: force=0, potential=0
                f_x(j) = 0;
                f_y(j) = 0;
                V_j(j) = 0;
            end
        end
        
        %Compute the position and velocity using Euler or Verlet
        if (strcmp('Kinematic',method)) %first timestap is always Euler   
            % Forward Euler x
            v(i,t+1) = sum(f_x)/lambda;
            x(i,t+1) = x(i,t) + v(i,t+1)*dt;
            
            % Forward Euler y
            u(i,t+1) = sum(f_y)/lambda;
            y(i,t+1) = y(i,t) + u(i,t+1)*dt;
            
        elseif(strcmp('Euler',method) || t==1)
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
        
        %%% If borders are on %%%%%%%%%%%%%%%       
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %Energy calculations
        V(i,t) = sum(V_j); %Potential energy: sum of influence of other particles
        T(i,t) = 1/2*m(i)*(v(i,t)^2+u(i,t)^2); %Kinetic energy: 1/2mv^2
    end
end




%%%%%%%%%%%%% Calculate and plot lumped quantities (energy calculated above) %%%%%%%%%%%%%%%%%%%%%%%%

for t = 1:timesteps 
    x_snapshot = x(:,t);
    y_snapshot = y(:,t);
    [k1,v1] = convhull(x(:,t),y(:,t));
    
    Z(:,t) = (1/N)*[sum(x(:,t)) sum(y(:,t))]'; % COM (\bar{x} \bar{y})
    X_bar = Z(1,:)';
    Y_bar = Z(2,:)';
    
    % fit ellipse and plot 
    ellipse_t = fit_ellipse(x_snapshot, y_snapshot);
        
    % grab ellipse parameters
    a = ellipse_t.long_axis;
    b = ellipse_t.short_axis;
    alpha = ellipse_t.phi;
    %orient = mod(alpha*180/pi,360);
    eccen = sqrt(1 - ((b/a)^2));
        
    % fill-in the ellipse values
    Alpha(t) = alpha;
    Orient(t) = alpha;
    Eccen(t) = eccen;
    A(t) = a;
    B(t) = b;
    Area(t) = v1;
     
end 



%%%%%%%%%%%%%%%%%%% Plot and save the figures %%%%%%%%%%%%%%%%%%%%%%%%%%%

% orientation
fig1 = figure(1);
plot(0:dt:time,Orient,'r','linewidth',3)
xlabel('Time','FontSize',14,'FontSize',14)
ylabel('Orientation','FontSize',14)
title('Orientation of the fit ellipse over time','FontSize',14)
fig1_name = 'Orientation.png';
saveas(fig1, fullfile(outpath, fig1_name))

% eccentricity
fig2 = figure(2);
plot(0:dt:time,Eccen,'k','linewidth',3)
xlabel('Time','FontSize',14)
ylabel('Eccentricity','FontSize',14)
title('Eccentricity of the ellipse varying over time','FontSize',14)
fig2_name = 'eccentricity.png';
saveas(fig2, fullfile(outpath, fig2_name))

% area
fig3 = figure(3);
plot(0:dt:time,Area,'b','linewidth',3)
xlabel('Time','FontSize',14)
ylabel('Area','FontSize',14)
title('Area of convex hull over time','FontSize',14)
fig3_name = 'Area.png';
saveas(fig3, fullfile(outpath, fig3_name))

% center of mass
fig4 = figure(4);
plot(X_bar,Y_bar,'k','linewidth',3)
title('Central of mass varying over time','FontSize',14)
xlabel('X','FontSize',14)
ylabel('Y','FontSize',14)
fig4_name = 'Center of Mass.png';
saveas(fig4, fullfile(outpath, fig4_name))

% energy
V_sum = sum(V)/2; %sum the potential energy of all particles times a half since matrix is mirrored
T_sum = sum(T) ;%1/2.*m'.*(sum(u.^2)+sum(v.^2)); %sum from all particles: 0.5*m*v^2
fig5=figure(5);
hold on
plot(0:dt:time,(T_sum+V_sum)/N^2)
plot(0:dt:time,T_sum/N^2,'r')
plot(0:dt:time,V_sum/N^2,'k')
title('The total energy of the system','FontSize',14);
xlabel('Time','FontSize',12);
ylabel('Energy','FontSize',12);
fig5_name = 'Energy.png';
saveas(fig5, fullfile(outpath, fig5_name))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%%%%%%%%%%%%%%%%%%% Draw the movie %%%%%%%%%%%%%%%%%%%%%

% set the axes 
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

% Draw the circles and make the movie
if ~(strcmp('', moviePath))
    fig=figure(6);
    mov = VideoWriter(moviePath, 'MPEG-4');

    for k=1:1/dt/fps:timesteps %pick the points for the correct number of fps
        circle = linspace(0,2*pi,100); %create a circle
        hold on;
        %plot N circles
        for n=1:N
            THETHA = atan2(u(n,round(k)),v(n,round(k)));
            xx = sigma/2*cos(circle)*cos(THETHA) - sigma/4*sin(circle)*sin(THETHA) + x(n,round(k)); %x-location of particle                  
            yy = sigma/2*cos(circle)*sin(THETHA) + sigma/4*sin(circle)*cos(THETHA) + y(n,round(k)); %y-location of particle

            if variables == 3
                color = [abs(v(n,round(k)))/max(max(abs(v))) 0 0];
            else
                color = [abs(v(n,round(k)))/max(max(abs(v))) 0 abs(u(n,round(k)))/max(max(abs(u)))];
            end
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
        THETHA_z = atan2(u_z(round(k)),v_z(round(k)));
        xx_z = sigma/2*cos(circle)*cos(THETHA_z) - sigma/4*sin(circle)*sin(THETHA_z) + x_z(round(k)); %x-location of particle                  
        yy_z = sigma/2*cos(circle)*sin(THETHA_z) + sigma/4*sin(circle)*cos(THETHA_z) + y_z(round(k)); %y-location of particle
        p3 = fill(xx_z,yy_z,'c');
        hold on
        
        % fit ellipse and plot 
        % NOTE: We have already fitted the ellipse above. Here we do it
        %       again to add it to the plot. Fix this
        ellipse_t = fit_ellipse(x_snapshot, y_snapshot, fig);
        
        % set the axes
        set(gca,'DataAspectRatio',[1 1 1]); %set aspect ratio x:y to 1:1
        axis([ax_x ax_y]); %set the axis
        %cleahold on;
        title('\fontsize{14} Circular motion of dog');
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
end




%%%%%%%%%%%%%%%%%% Save the workspace %%%%%%%%%%%%%%%%%%%%%%%%%%%
workspacePath = fullfile(outpath, 'Variables');
save(workspacePath, 'sigma', 'epsilon', 'L', 'l', 'G', 'lambda', 'c', 'time', 'dt', 'method', 'borders', 'ApplyBC', 'fps', ...
                    'N', 'm', 'timesteps', 'V', 'T', 'x', 'y', 'u', 'v', 'f_x', 'f_y', 'V_j', 'A', 'B', 'save_Z', 'save_X_bar', ...
                    'save_Y_bar', 'Alpha', 'Orient', 'Eccen', 'Z', 'X_bar', 'Y_bar', 'Area');


