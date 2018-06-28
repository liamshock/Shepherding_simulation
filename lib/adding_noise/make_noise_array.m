function [noise_array] = make_noise_array(numSheep, initCond, bias, stdev, memoryTime, dt, finalTime)
% Make the noise array
% 
% -- Output --
% noise_array: shape = (timesteps, numSheep, 2), where last dim is X and Y vals

% get the number of timesteps
timesteps = round(finalTime/dt)+1;

% allocate the output array
noise_array = zeros(timesteps, numSheep, 2);

% fill in the noise array
for i = 1:numSheep
    noise_array(:, i, 1) = pinknoise(initCond, bias, stdev, memoryTime, dt, finalTime);
    noise_array(:, i, 2) = pinknoise(initCond, bias, stdev, memoryTime, dt, finalTime);
end

end

