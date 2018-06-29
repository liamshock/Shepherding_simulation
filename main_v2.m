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
c = 1*ones(N_dogs,1);                   % predator-prey repulsion
maxv = 2;                                 % max velocity of sheep
time = 40 ;                                 % total time
dt = 0.001;                              % time step
method = 'Kinematic';                     % integration method
borders = [5 10 20 25];                      % set [x0 x1 y0 y1] to turn borders on (particles bounce), 0 is off
ApplyBC = false;                          % should we use the boundary conditions
N = 50;                                   % num of particles 
fps = 10;                                 % FPS for movie
% x_T = 10;                                 % desired x-comp of flock COM
% y_T = 10;                                 % desired x-comp of flock COM
spd_dog = 4*ones(N_dogs,1);               % speed of the dog
% beta = pi/3;                              % dog driving angle
% dog_dist = sqrt(N).*ones(N_dogs,1);       % the distance the dog keeps from the herd
tau = 0.01/dt;                               % num of timesteps we wait for flock relaxation
DBeta = pi/6;
tol = 0.4;
stdev_noise = 2.5;
memoryTime_noise = 1;
dog_rad_fac = 1.5;

% % load the trajectories
load Targets.mat X_T Y_T 

% garden 
Wx = 50; Wy = 30; d1 = 10; d2 = 12; d3 = 2; d4 = 6; 




% % targeted motion 
% Wx = 50; Wy = 30; d1 = 10; d2 = 12; d3 = 2; d4 = 6; 
% Fullborders = [0 Wx 0 Wy];                      
% gpx1 = Wx-d2; gpx2 = gpx1+0.5*(d2-d4); gpx3 = gpx2+d4; gpx4 = Wx; gpx5 = gpx1; gpx6 = gpx1-d3;
% gpy1 = Wy-d1; gpy2 = gpy1; gpy3 = gpy1; gpy4 = gpy1; gpy5 = Wy; gpy6 = gpy5;
% x_T = 0.5*(gpx1+gpx2);                                 % desired x-comp of flock COM
% y_T = 0.5*(gpy1+gpy6);                                % desired x-comp of flock COM
% 
% % inspection figure 
% figure 
% plot([gpx1,gpx6],[gpy1,gpy6],'k','linewidth',3);hold on
% plot([gpx1,gpx2],[gpy1,gpy2],'k','linewidth',3);
% plot([gpx3,gpx4],[gpy3,gpy4],'k','linewidth',3);
% plot([0,Wx],[0,0],'k','linewidth',3);
% plot([Wx,Wx],[0,Wy],'k','linewidth',3);
% plot([Wx,0],[Wy,Wy],'k','linewidth',3);
% plot([0,0],[Wy,0],'k','linewidth',3);
% axis equal tight 
% hold on
% plot([borders(1:2)],[borders(3),borders(3)],'r','linewidth',2);
% plot([borders(2),borders(2)],[borders(3:4)],'r','linewidth',2);
% plot([borders(2),borders(1)],[borders(4),borders(4)],'r','linewidth',2);
% plot([borders(1),borders(1)],[borders(4),borders(3)],'r','linewidth',2);
% 
% xpthIc = mean(borders(1:2));
% ypthIc = mean(borders(3:4));
% 
% 
% plot(xpthIc,ypthIc,'*r');
% plot(x_T,y_T,'*r');
% % [x,y] = ginput(7);
% load TargTrajPts.mat x y 
% 
% plot(x,y,'--r','linewidth',2);







%% %%%%%%%%%%%%%%%% Main body %%%%%%%%%%%%%%%%%%%%%%%%%

% run the main simulation
[x, y, u, v, f_x, f_y, V_j, V, T, x_dog, y_dog, u_dog, v_dog, noise_arr] = run_simulation(maxv, tau, X_T, Y_T, spd_dog, ...
                                                                                                                  seed, epsilon, c, time, dt, N, mass, ...
                                                                                                                  method, borders, ApplyBC, ...
                                                                                                                  N_dogs, DBeta, tol, ...
                                                                                                                  Wx, Wy, d1, d2, d3, d4, ...
                                                                                                                  stdev_noise, memoryTime_noise, ...
                                                                                                                  dog_rad_fac);
                                          
% extract the data arrays
[A, B, orientation, eccentricity, X_bar, Y_bar, area] = extract_data(x, y, time, dt, N);

% plot the functions
plot_functions(outpath, time, dt, N, orientation, eccentricity, area, X_bar, Y_bar, T, V);

% write the movie
make_movie(outpath, dt, time, fps, N, ApplyBC, borders, x, y, u, v, x_dog, y_dog, u_dog, v_dog, X_T, Y_T, N_dogs, DBeta, tol,Wx,Wy, d1, d2, d3, d4);

% save the workspace 
workspacePath = fullfile(outpath, 'Variables');
save(workspacePath, 'epsilon', 'c', 'time', 'dt', 'method', 'borders', 'ApplyBC', 'fps', ...
                    'N', 'mass', 'V', 'T', 'x', 'y', 'u', 'v', 'f_x', 'f_y', 'V_j', 'A', 'B', ...
                    'orientation', 'eccentricity', 'X_bar', 'Y_bar', 'area');
                


