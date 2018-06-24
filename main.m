% This is the main run script for the simulation using functions


% -- Main functions called -- %
% run_simulation -> evolve the position of the sheep
% extract_data   -> extract the reduced parameters from the simulation
%                   results
% plot_functions -> save plots of the reduced parameters
% make_movie     -> make and save a movie of the simulation


% -- Notes -- %
% (1): For now we prescribe the path of the dog and pass it directly to the
%      simulation
% (2): We are supressing the output of the images and movies. They are
%      saved directly



%%%%%%%%%%%%%%%%%   Parameters %%%%%%%%%%%%%%%%%%%%%%

outpath = '/Users/liam/output';           % the directory where information will be saved
seed = 4;                                 % the seed for the random number generator
sigma = 1;                                % particle diameter
epsilon = 15;                             % material parameter
L = 20;                                   % attractive length in particle diameters
l = 5;                                    % repulsive length from the predator
G = 0.1;                                  % attractive strength
lambda = 1;                               % damping coefficient
c = 0.4;                                  % predator-prey repulsion
time = 2;                                 % total time
dt = 0.001;                               % time step
method = 'Kinematic';                     % integration method
borders = [-30 30 -30 30];                % set [x0 x1 y0 y1] to turn borders on (particles bounce), 0 is off
ApplyBC = false;                          % should we use the boundary conditions
N = 20;                                   % num of particles 
mass = 0.1;                               % 'social mass' of the sheep
fps = 10;                                 % FPS for movie




%%%%%%%%%%%%%%% Prescribe the dog motions %%%%%%%%%%%%%

% Init Position and velocity of a single Predator
thet0_z = linspace(0,4*pi,time/dt+1); 
x_dog = 5.*cos(thet0_z)+15; 
y_dog = 5.*sin(thet0_z)+15;
v_dog = -5.*sin(thet0_z); 
u_dog = 5.*cos(thet0_z);




%%%%%%%%%%%%%%%%%% Main body %%%%%%%%%%%%%%%%%%%%%%%%%

% run the main simulation
[x, y, u, v, f_x, f_y, V_j, T, V] = run_simulation(seed, sigma, epsilon, L, l, G, ...
                                                   lambda, c, time, dt, N, mass, ...
                                                   method, borders, ApplyBC, ...
                                                   x_dog, y_dog);
                                          
% extract the data arrays
[A, B, orientation, eccentricity, X_bar, Y_bar, area] = extract_data(x, y, time, dt, N);

% plot the functions
plot_functions(outpath, time, dt, N, orientation, eccentricity, area, X_bar, Y_bar, T, V)

% write the movie
make_movie(outpath, sigma, dt, time, fps, N, ApplyBC, borders, x, y, u, v, x_dog, y_dog, u_dog, v_dog)

% save the workspace 
workspacePath = fullfile(outpath, 'Variables');
save(workspacePath, 'sigma', 'epsilon', 'L', 'l', 'G', 'lambda', 'c', 'time', 'dt', 'method', 'borders', 'ApplyBC', 'fps', ...
                    'N', 'mass', 'V', 'T', 'x', 'y', 'u', 'v', 'f_x', 'f_y', 'V_j', 'A', 'B', ...
                    'orientation', 'eccentricity', 'X_bar', 'Y_bar', 'area');


