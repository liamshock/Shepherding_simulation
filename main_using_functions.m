% This is the main run script for the simulation using functions


% -- Main functions called -- %
% run_simulation -> evolve the position of the sheep


% -- Notes -- %
% (1): For now we 

%%%%%%%%%%%%%%%%%   Parameters %%%%%%%%%%%%%%%%%%%%%%

seed = 4;                  % the seed for the random number generator
sigma = 1;                 % particle diameter
epsilon = 15;              % material parameter
L = 20;                    % attractive length in particle diameters
l = 5;                     % repulsive length from the predator
G = 0.1;                   % attractive strength
lambda = 1;                % damping coefficient
c = 0.4;                   % predator-prey repulsion
time = 1;                 % total time
dt = 0.001;                % time step
method = 'Kinematic';      % integration method
borders = [-30 30 -30 30]; % set [x0 x1 y0 y1] to turn borders on (particles bounce), 0 is off
ApplyBC = false;           % should we use the boundary conditions
N = 50;                    % num of particles 
mass = 0.1;                % mass of the particles




%%%%%%%%%%%%%%% Prescribe the dog motions %%%%%%%%%%%%%

% Init Position and velocity of a single Predator
thet0_z = linspace(0,4*pi,time/dt+1); 
x_dog = 5.*cos(thet0_z)+15; 
y_dog = 5.*sin(thet0_z)+15;
v_dog = -5.*sin(thet0_z); 
u_dog = 5.*cos(thet0_z);




%%%%%%%%%%%%%%%%%% Main body %%%%%%%%%%%%%%%%%%%%%%%%%

% run the main simulation
[x, y, u, v, f_x, f_y, V_j, T, V] = run_calculation(seed, sigma, epsilon, L, l, G, ...
                                              lambda, c, time, dt, N, mass, ...
                                              method, borders, ApplyBC, ...
                                              x_dog, y_dog);
                                          
% extract the data arrays
[A, B, Orient, Eccen, X_bar, Y_bar, Area] = extract_data(x, y);


