function [x_dog_new, y_dog_new] = dog_position_v2(xd_curr,yd_curr,xbar, ybar, phi, N_dogs, dog_dist, x_T, y_T, maxDis, DBeta)
% find the position of the dog that would move the herd in the correct
% direction

    Beta = atan2((-ybar+y_T),(-xbar+x_T));
    Edel = (2*pi-(DBeta*DBeta+Beta))/4;
    Th1 = Beta+DBeta+Edel;
    Th2 = Th1+Edel;
    Th3 = Th2+Edel;
    radDogDist = 1.1*maxDis;
    px1 = xbar + radDogDist*cos(Th1); py1 = ybar + radDogDist*sin(Th1);
    px2 = xbar + radDogDist*cos(Th2); py2 = ybar + radDogDist*sin(Th2);
    px3 = xbar + radDogDist*cos(Th3); py3 = ybar + radDogDist*sin(Th3);
    targetXdog = [px1,px2,px3]; targetYdog = [py1,py2,py3];
    
    DDMat = zeros(N_dogs,N_dogs);
    DogDestIdx = nan(N_dogs,1);
    for kkdog = 1:N_dogs
    DDMat(kkdog,:) = sqrt((xd_curr(kkdog)-targetXdog).^2+(yd_curr(kkdog)-targetYdog).^2);
    end
        
    for kkdog = 1:N_dogs
    [~,idx] = min(DDMat(kkdog,:));
    DogDestIdx(kkdog) = idx;
    DDMat(:,idx) = 1000;
    end
    
    x_dog_new = targetXdog(DogDestIdx); x_dog_new=x_dog_new';
    y_dog_new = targetYdog(DogDestIdx); y_dog_new=y_dog_new';

end

