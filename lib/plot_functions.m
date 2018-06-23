function [] = plot_functions(outpath, time, dt, N, orientation, eccentricity, area, X_bar, Y_bar, Kinetic_E, Potential_E)
% Save plots to the output path


% Introductory message
fprintf('\n\n')
fprintf('--- Making and saving the plots ---')  
fprintf('\n')

% First we create the output directory if it does not already exist
if ~exist(outpath, 'dir')
  mkdir(outpath);
end

% orientation
fig = figure('visible', 'off');
plot(0:dt:time, orientation,'r','linewidth',3)
xlabel('Time','FontSize',14,'FontSize',14)
ylabel('Orientation','FontSize',14)
title('Orientation of the fit ellipse over time','FontSize',14)
fig_name = 'Orientation.png';
saveas(fig, fullfile(outpath, fig_name))
close(fig)

% eccentricity
fig = figure('visible', 'off');
plot(0:dt:time, eccentricity,'k','linewidth',3)
xlabel('Time','FontSize',14)
ylabel('Eccentricity','FontSize',14)
title('Eccentricity of the ellipse varying over time','FontSize',14)
fig_name = 'eccentricity.png';
saveas(fig, fullfile(outpath, fig_name))
close(fig)

% area
fig = figure('visible', 'off');
plot(0:dt:time, area,'b','linewidth',3)
xlabel('Time','FontSize',14)
ylabel('Area','FontSize',14)
title('Area of convex hull over time','FontSize',14)
fig_name = 'Area.png';
saveas(fig, fullfile(outpath, fig_name))
close(fig)

% center of mass
fig = figure('visible', 'off');
plot(X_bar,Y_bar,'k','linewidth',3)
title('Central of mass varying over time','FontSize',14)
xlabel('X','FontSize',14)
ylabel('Y','FontSize',14)
fig_name = 'Center of Mass.png';
saveas(fig, fullfile(outpath, fig_name))
close(fig)

% energy
V_sum = sum(Potential_E)/2; %sum the potential energy of all particles times a half since matrix is mirrored
T_sum = sum(Kinetic_E) ;%1/2.*m'.*(sum(u.^2)+sum(v.^2)); %sum from all particles: 0.5*m*v^2
fig=figure('visible', 'off');
hold on
plot(0:dt:time,(T_sum+V_sum)/N^2)
plot(0:dt:time,T_sum/N^2,'r')
plot(0:dt:time,V_sum/N^2,'k')
title('The total energy of the system','FontSize',14);
xlabel('Time','FontSize',12);
ylabel('Energy','FontSize',12);
fig_name = 'Energy.png';
saveas(fig, fullfile(outpath, fig_name))
close(fig)

% end of the function
fprintf('Plots saved')
fprintf('\n') 
end

