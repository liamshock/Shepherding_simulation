% This is the main run script for the simulation using functions

clear 
close all
clc

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

% the directory where information will be saved
outpath = strcat('/Users/liam/output/', datestr(datetime('now')));   

% set the directory where the code is being stored (we will save a copy)
codeDirectory = '/Users/liam/code/Shepherding_simulation';

% save a copy of the code as it currently stands
save_all_code(codeDirectory, outpath);


seed = 4;                                 % the seed for the random number generator
N_dogs = 3;                               % num of dogs
epsilon = 1;                              % material parameter
mass = 0.1;                               % 'social mass' of the sheep
c = 0.0*ones(N_dogs,1);                   % predator-prey repulsion
maxv = 1;                                 % max velocity of sheep
time = 10;                                 % total time
dt = 0.0001;                              % time step
method = 'Kinematic';                     % integration method
borders = [0 5 0 5];                    % set [x0 x1 y0 y1] to turn borders on (particles bounce), 0 is off
ApplyBC = false;                          % should we use the boundary conditions
N = 30;                                   % num of particles 
fps = 10;                                 % FPS for movie
x_T = 10;                                 % desired x-comp of flock COM
y_T = 10;                                 % desired x-comp of flock COM
spd_dog = 2*ones(N_dogs,1);               % speed of the dog
beta = pi/3;                              % dog driving angle
dog_dist = sqrt(N).*ones(N_dogs,1);       % the distance the dog keeps from the herd
tau = 9.9/dt;                               % the number of timesteps we wait for equilibrium of flock (should be less than time)





%%%%%%%%%%%%%%%%%% Main body %%%%%%%%%%%%%%%%%%%%%%%%%

% run the main simulation
[x, y, u, v, f_x, f_y, V_j, V, T, x_dog, y_dog, u_dog, v_dog, x_bar_init, y_bar_init, noise_arr] = run_simulation(maxv, tau, x_T, y_T, spd_dog, beta, seed, ...
                                                                                                       epsilon, c, time, dt, N, mass, ...
                                                                                                       method, borders, ApplyBC, dog_dist, N_dogs);
                                          
% extract the data arrays
[A, B, orientation, eccentricity, X_bar, Y_bar, area] = extract_data(x, y, time, dt, N);

% plot the functions
plot_functions(outpath, time, dt, N, orientation, eccentricity, area, X_bar, Y_bar, T, V);

% write the movie
make_movie(outpath, dt, time, fps, N, ApplyBC, borders, x, y, u, v, x_dog, y_dog, u_dog, v_dog, x_bar_init, y_bar_init, x_T, y_T, N_dogs);

% save the workspace 
workspacePath = fullfile(outpath, 'Variables');
save(workspacePath, 'epsilon', 'c', 'time', 'dt', 'method', 'borders', 'ApplyBC', 'fps', ...
                    'N', 'mass', 'V', 'T', 'x', 'y', 'u', 'v', 'f_x', 'f_y', 'V_j', 'A', 'B', ...
                    'orientation', 'eccentricity', 'X_bar', 'Y_bar', 'area');
                


