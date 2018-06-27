function noise = pinknoise(initCond, bias, stdev, memoryTime, dt, finalTime, seed, figNumber)

    if exist('seed','var')
        rng(seed)
    end
    t = dt:dt:finalTime;      % Time vector
    
    ex = exp(-t/memoryTime);
    noise = initCond*ex+bias*(1-ex)+stdev*ex.*cumsum([0 sqrt(diff(1./ex.^2-1)).*randn(1,length(t)-1)])/sqrt(2/memoryTime);
    
    if exist('figNumber','var')
        figure(figNumber);
        plot(t,noise);
    end

end