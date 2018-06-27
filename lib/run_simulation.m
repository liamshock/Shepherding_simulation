function [x, y, u, v, f_x, f_y, V_j, V, T, x_dog, y_dog, u_dog, v_dog, x_bar_init, y_bar_init] = run_simulation(maxv, tau, x_T, y_T, spd_dog, beta, seed, sigma, epsilon, L, l, G, ...
                                                                                                               lambda, c, time, dt, N, mass, ...
                                                                                                               method, borders, ApplyBC, dog_dist, N_dogs)
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
m = mass*ones(N,1);       %used for Euler or Verlet only

%number of timesteps to be calculated
timesteps = round(time/dt)+1;


% initialize variables 

% sheep
x = zeros(N,timesteps); 
y = x;  
u = x; 
v = x;

% Initial random positions for sheep
x0 = borders(1) + (borders(2)-borders(1))*rand(N,1); 
y0 = borders(3) + (borders(4)-borders(3))*rand(N,1);

% dogs
x_dog = zeros(N_dogs, timesteps);
y_dog = x_dog;
v_dog = x_dog;
u_dog = x_dog;

% Initial random positions for dogs (first frame)
% NOTE: FOR ONE DOG
% x_dog(:, 1) = -10*ones(N_dogs,1);
% y_dog(:, 1) = -10*ones(N_dogs,1);
x_dog(:, 1) = [-2;-1;0];
y_dog(:, 1) = [-2;-1;0];
v_dog(:, 1) = zeros(N_dogs,1);
u_dog(:, 1) = zeros(N_dogs,1);


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

riz = zeros(N_dogs,1);
n_dog = zeros(N_dogs,2);
% f_xdog = zeros(1,N);
% f_ydog = zeros(1,N);
V_jdog = zeros(1,N);


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
                for k = 1:N_dogs
                    riz(k) = sqrt((x(i,t) - x_dog(k,t))^2+(y(i,t)-y_dog(k,t))^2);
                end
                
                %compute unit normal vector, n(1)=n_x, n(2)=n_y
                n = [ x(j,t)-x(i,t) y(j,t)-y(i,t) ]/norm([ x(j,t) - x(i,t) y(j,t)- y(i,t) ]);
                
                for k = 1:N_dogs
                    n_dog (k,:)= [ x_dog(k,t)-x(i,t) y_dog(k,t)-y(i,t)]/norm([ x_dog(k,t) - x(i,t) y_dog(k,t) - y(i,t)]);
                end
                
                f_xdog_timestep = zeros(N_dogs, 1);
                f_ydog_timestep = zeros(N_dogs, 1);
                V_jdog_timestep = zeros(N_dogs, 1);
                for k = 1:N_dogs
%                     for x direction: force between Prey & predator
                    %f_xdog_timestep(k) = (n_dog(k,1)*c(k)*exp(-riz(k)/(l(k)*sigma))/(l(k)*sigma));
                    f_xdog_timestep(k) = n_dog(k,1)*c(k)/riz(k);

%                     for y direction: force between Prey & predator
                    %f_ydog_timestep(k) = (n_dog(k,2)*c(k)*exp(-riz(k)/(l(k)*sigma))/(l(k)*sigma));
                    f_ydog_timestep(k) = n_dog(k,2)*c(k)/riz(k); 
                    
%                     potential energy for dogs only
                    %V_jdog(k) =  V_jdog(k) - c(k)*exp(-riz(k)/(l(k)*sigma));
                    V_jdog_timestep(k) = c(k)*log(abs(riz(k)));
%                 
                end
                dog_forces_x = sum(f_xdog_timestep);
                dog_forces_y = sum(f_ydog_timestep);
                dog_potentials = sum(V_jdog_timestep);
                
                %for x direction
                %f_x(j) = -4*epsilon*(n(1)*exp(-rij/sigma)/sigma - n(1)*G*exp(-rij/(L*sigma))/(L*sigma) + dog_forces_x);
                f_x(j) = -((n(1)/rij - n(1)*rij) + dog_forces_x);
                
                %for y direction
                %f_y(j) = -4*epsilon*(n(2)*exp(-rij/sigma)/sigma - n(2)*G*exp(-rij/(L*sigma))/(L*sigma) + dog_forces_y);
                f_y(j) = -((n(2)/rij - n(2)*rij) + dog_forces_y);
                
                %potential energy
                %V_j(j) = 4*epsilon*(-sigma*exp(-rij/sigma) + G*exp(-rij/(L*sigma)) - c(1)*exp(-riz(1)/(l(1)*sigma)) + V_jdog(j));
                V_j(j) = (log(abs(rij)) - 0.5*rij^2) + dog_potentials;
                
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
            x(i,t+1) = x(i,t) + sign(v(i,t+1))*min(abs(v(i,t+1)), maxv)*dt;
            
            % Forward Euler y
            u(i,t+1) = sum(f_y)/lambda;
            y(i,t+1) = y(i,t) + sign(u(i,t+1))*min(abs(u(i,t+1)), maxv)*dt;
            
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
    
    
    % define theta after tau steps
    if t == tau
        z_bar = (1/N)*[sum(x(:,t)) sum(y(:,t))]'; % COM (\bar{x} \bar{y})
        x_bar_init = z_bar(1);
        y_bar_init = z_bar(2);
        theta = atan2(y_T - y_bar_init, x_T - x_bar_init);
    end
    

    
    % after transients
    if t > tau
        
       
        % above or below the line? Use cross product
        x_bar = (1/N)*(sum(x(:,t)));
        y_bar = (1/N)*sum(y(:,t));
        z_bar = [x_bar y_bar];   
        CM_T = [x_T - x_bar_init, y_T - y_bar_init];
        stack = [z_bar; CM_T];
        sign_det = sign(det(stack));
    
        % define phi
        if sign_det < 0
            phi = (3*pi/2) + theta + beta;
        else
            phi = theta + (pi/2) - beta;
        end
        
        % calculate error
        %error_modulus = find_error_modulus(x_bar_init, y_bar_init, x_bar, y_bar, theta);
      
        % find the desired position of the dog
        [x_dog_desired, y_dog_desired] = dog_position(x_bar, y_bar, phi, N_dogs, dog_dist);
                      
%         % update dog position and velocity
        for k = 1:N_dogs
            v_dog(k, t+1) = ((x_dog_desired(k) - x_dog(k,t)) / norm(x_dog_desired(k) - x_dog(k,t)))*spd_dog(k);
            u_dog(k, t+1) = ((y_dog_desired(k) - y_dog(k,t)) / norm(y_dog_desired(k) - y_dog(k,t)))*spd_dog(k);
            x_dog(k, t+1) = x_dog(k, t) + v_dog(k, t+1)*dt;
            y_dog(k, t+1) = y_dog(k, t) + u_dog(k, t+1)*dt;
        end
    end
    
end


% end of the function
fprintf('Simulation has finished running')
fprintf('\n')
end

