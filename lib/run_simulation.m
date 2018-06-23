function [x, y, u, v, f_x, f_y, V_j, V, T] = run_simulation(seed, sigma, epsilon, L, l, G, ...
                                                       lambda, c, time, dt, N, mass, ...
                                                       method, borders, ApplyBC, ...
                                                       x_dog, y_dog)
% Run the main simulation
% We evolve the position of the sheep over time by applying the forces
% exerted on each other and the forces exerted by the dog
                                                   

% introductory message
fprintf('\n\n')
fprintf('--- Running the simulation ---')  
fprintf('\n')

% set the seed
rng(seed)
                
% assign the mass of the particles
m = mass*ones(N,1); 

%number of timesteps to be calculated
timesteps = time/dt+1;

% initialize variables 
x = zeros(N,timesteps); 
y = x;  
u = x; 
v = x;

% Initial random positions for sheep
x0 = borders(1) + (borders(2)-borders(1))*rand(N,1); 
y0 = borders(3) + (borders(4)-borders(3))*rand(N,1);

% initial velocities 
thet0 = 2*pi*rand(N,1); 
spd0 = rand(N,1); 
u0 = spd0.*cos(thet0);
v0 = spd0.*sin(thet0);

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


%%%%%%%%%%%%%%%%%%% Calcualtion part %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for t=1:timesteps
    
        % print quarterly updates
    if t == round(timesteps/10)
        fprintf('1/10 complete \n')
    elseif t == round(timesteps/5)
        fprintf('2/10 complete \n')
    elseif t == round( (3/10)*timesteps)
        fprintf('3/10 complete \n')
    elseif t == round( (4/10)*timesteps)
        fprintf('4/10 complete \n')   
    elseif t == round(timesteps/2)
        fprintf('5/10 complete \n')   
    elseif t == round( (6/10)*timesteps)
        fprintf('6/10 complete \n')  
    elseif t == round( (7/10)*timesteps)
        fprintf('7/10 complete \n')
    elseif t == round( (8/10)*timesteps)
        fprintf('8/10 complete \n')
    elseif t == round( (9/10)*timesteps)
        fprintf('9/10 complete \n')
    end
    
    for i=1:N
        for j=1:N %set force/energy interactions for each particle
            if i~=j %don't take self inducting terms into account
                
                %compute the distance
                rij = sqrt((x(i,t) - x(j,t))^2+ (y(i,t)-y(j,t))^2);
                riz = sqrt((x(i,t) - x_dog(t))^2+(y(i,t)-y_dog(t))^2);
                
                %compute unit normal vector, n(1)=n_x, n(2)=n_y
                n = [ x(j,t)-x(i,t) y(j,t)-y(i,t) ]/norm([ x(j,t) - x(i,t) y(j,t)- y(i,t) ]);
                n_z = [ x_dog(t)-x(i,t) y_dog(t)-y(i,t) ]/norm([ x_dog(t) - x(i,t) y_dog(t) - y(i,t) ]);
                
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


% end of the function
fprintf('Simulation has finished running')
fprintf('\n')
end

