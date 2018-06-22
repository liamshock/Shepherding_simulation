%This script plots a bunch of particles and their behavior
%The script can read both 1D and 2D files

close all
clear all
clc

%%% DATA TO CHANGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rng(4);
sigma = 1;              %particle diameter
epsilon = 15;           %material parameter
L = 30;                 %attractive length in particle diameters
l = 10;                 %repulsive length from the predator
G = 0.5;                  %attractive strength
lambda = 1;           %damping coefficient
c = 0.27;                % predator-prey repulsion
%rep = -1;               %repulsion from shepherd (set -1 for no shepherd)
time = 50;              %total time
dt = 0.001;             %time step
method = 'Kinematic';      %integration method
borders = [00 30 00 30];%set [x0 x1 y0 y1] to turn borders on (particles bounce), 0 is off
output = '';            %name of the output file. Set '' for no output
movie = 'Ellip_circl_pred_fitE.mp4';   %movie output name. Set '' for no movie output
fps = 10;               %FPS for movie
ApplyBC = false;

N = 30; % num of particles 
m = 0.1*ones(N,1); % mass of the particles

timesteps = time/dt+1; %number of timesteps to be calculated

% initialize variables 
x = zeros(N,timesteps); y = x;  u = x; v = x;

%%%%%%%%%%%%%
thet0_z = linspace(0,4*pi,timesteps); %2*pi*rand(1);
%spd0_z = 2;             % Speed of the Predator
% u0_z = cos(thet0_z);
% v0_z = sin(thet0_z);

% x_z = zeros(1,timesteps)+15; y_z = x_z;    % Init Position of a single Predator
x_z = 5.*cos(thet0_z)+15; y_z = 5.*sin(thet0_z)+15;
% x_z = 15 ; y_z=15;
v_z = -5.*sin(thet0_z); u_z = 5.*cos(thet0_z);

% Initial conditions and velocities 
x0 = borders(1) + borders(2)*rand(N,1); % x - initial positions of agents form uniform random distr 
y0 = borders(3) + borders(4)*rand(N,1); % y - initial positions of agents form uniform random distr 

% initial velocities 
thet0 = 2*pi*rand(N,1); 
spd0 = rand(N,1); 
u0 = spd0.*cos(thet0);
v0 = spd0.*sin(thet0);
variables = 5;

% figure 
% plot(x0,y0,'.k')

% place the Ic in their corresponding matrices 
x(:,1) = x0;
y(:,1) = y0;
u(:,1) = u0;
v(:,1) = v0;
% x_z(1) = x0_z;
% y_z(1) = y0_z;

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
Area_ell = zeros([timesteps 1]);
%Area_conv = zeros([timesteps 1]);




%%% Calcualtion part %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t=1:time/dt+1
    for i=1:N
        for j=1:N %set force/energy interactions for each particle
            if i~=j %don't take self inducting terms into account
                
                %compute the distance
                rij = sqrt((x(i,t)-x(j,t))^2+ (y(i,t)-y(j,t))^2);
                riz = sqrt((x(i,t)-x_z(t))^2+(y(i,t)-y_z(t))^2);
                
                %compute unit normal vector, n(1)=n_x, n(2)=n_y
                n = [ x(j,t)-x(i,t) y(j,t)-y(i,t) ]/norm([ x(j,t)-x(i,t) y(j,t)-y(i,t) ]);
                
                %compute effective attractive strength
                %Geff = G*(1+(j==1)*(-rep-1));
                
                %for x direction
               %f_x(j) = n(1)*-4*epsilon*((6*R*(L*sigma)^6)/rij^7 - (12*(L*sigma)^12)/rij^13);
                f_x(j) = n(1)*-4*epsilon*(exp(-rij/sigma)/sigma - G*exp(-rij/(L*sigma))/(L*sigma)+c*exp(-riz/(l*sigma))/(l*sigma));
                %f_x(j) = n(1)*-4*epsilon*(1/rij) - Geff*(rij);
                
                %for y direction
               %f_y(j) = n(2)*-4*epsilon*((6*(L*sigma)^6)/rij^7 - (12*(L*sigma)^12)/rij^13);
                f_y(j) = n(2)*-4*epsilon*(exp(-rij/sigma)/sigma - G*exp(-rij/(L*sigma))/(L*sigma)+c*exp(-riz/(l*sigma))/(l*sigma));
                %f_y(j) = n(2)*-4*epsilon*(1/rij) - Geff*(rij);
                
                %potential energy
               %V_j(j) = 4*epsilon*((sigma/rij)^12-(L*sigma/rij)^6);
                V_j(j) = 4*epsilon*(G*exp(-rij/(L*sigma)) + c*exp(-riz/(l*sigma)) - exp(-rij/sigma));
                %V_j(j) = 4*epsilon*(log(abs(rij)) - Geff*(1/2)*((rij)^2));
            
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
    Z(:,t) = (1/N)*[sum(x(:,t)) sum(y(:,t))]'; % COM (\bar{x} \bar{y})
    X_bar = Z(1,:)';
    Y_bar = Z(2,:)';
end

% Calculate and plot lumped quantities %%%%%%%%%%%%%%%%%%%%%%%%

for t = 1:timesteps 
    pos_curr = [x(:,t) y(:,t)]';
    % Fit ellipse to data 
    [z, a, b, alpha] = fitellipse(pos_curr);
    Q = [cos(alpha), -sin(alpha); 
         sin(alpha), cos(alpha)];
    theta = linspace(0,2*pi,360);
    PosEllip = z + Q*[a * cos(theta); b * sin(theta)];
    orient = mod(alpha*180/pi,360);
    eccen = b/a;
    
    A(t) = a;
    B(t) = b;
    save_Z(t,:) = z;
    save_X_bar(t) = z(1);
    save_Y_bar(t) = z(2);
    Alpha(t) = alpha;
    Orient(t) = orient;
    Eccen(t) = eccen;
    Area_ell(t) = pi*a*b;
    
end 


figure(1)
plot(0:dt:time,Orient,'r','linewidth',3)
xlabel('Time','FontSize',14,'FontSize',14)
ylabel('Orientation','FontSize',14)
title('Orientation of the fit ellipse over time','FontSize',14)
figure(2)
plot(0:dt:time,Eccen,'k','linewidth',3)
xlabel('Time','FontSize',14)
ylabel('Eccentricity','FontSize',14)
title('Eccentricity of the ellipse varying over time','FontSize',14)
figure(3)
plot(0:dt:time,Area_ell,'b','linewidth',3)
xlabel('Time','FontSize',14)
ylabel('Area','FontSize',14)
title('Area using Ellipse fit over time','FontSize',14)
figure(4)
plot(X_bar,Y_bar,'k','linewidth',3)
title('Central of mass varying over time','FontSize',14)
xlabel('X','FontSize',14)
ylabel('Y','FontSize',14)


% Calculate and plot energy %%%%%%%%%%%%%%%%%%%%%%%%
V_sum = sum(V)/2; %sum the potential energy of all particles times a half since matrix is mirrored
T_sum = sum(T) ;%1/2.*m'.*(sum(u.^2)+sum(v.^2)); %sum from all particles: 0.5*m*v^2

fig=figure(5);
hold on
plot(0:dt:time,(T_sum+V_sum)/N^2)
plot(0:dt:time,T_sum/N^2,'r')
plot(0:dt:time,V_sum/N^2,'k')
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
    fig=figure(6);
    mov = VideoWriter(movie, 'MPEG-4');

    for k=1:1/dt/fps:timesteps %pick the points for the correct number of fps
        circle = linspace(0,2*pi,100); %create a circle
        hold on;
        for n=1:N %plot N circles...
%             xx = sigma/2*cos(circle)+x(n,round(k)); %x-location of particle                  
%             yy = sigma/2*sin(circle)+y(n,round(k)); %y-location of particle
            THETHA = atan2(u(n,round(k)),v(n,round(k)));
            %THETHA = 1
%             dott = dot(u(n,round(k)),v(n,round(k)));
%             crosss = norm(cross(u(n,round(k)),v(n,round(k))));
%             THETHA = atan(crosss/dott);
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
%         xx_z = sigma/2*cos(circle)+x_z(round(k)); %x-location of particle                  
%         yy_z = sigma/2*sin(circle)+y_z(round(k)); %y-location of particle
        p3 = fill(xx_z,yy_z,'c');
        hold on
        
        % fit ellipse
        plotellipse(save_Z(round(k),:), A(round(k)), B(round(k)), Alpha(round(k)))
      
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
    %mov = close(mov);
end
close(mov);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




