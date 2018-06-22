% fitting the ellipse

x_snap = x(:,1)'
y_snap = y(:,1)'

fig1 = figure(1)


ellipse_t = fit_ellipse(x_snap, y_snap, fig1)

fig1 = figure(1)


        circle = linspace(0,2*pi,100); %create a circle
        hold on;
        %plot N circles
        for n=1:20
            THETHA = atan2(u(n,round(k)),v(n,round(k)));
            xx = sigma/2*cos(circle)*cos(THETHA) - sigma/4*sin(circle)*sin(THETHA) + x(n,round(k)); %x-location of particle                  
            yy = sigma/2*cos(circle)*sin(THETHA) + sigma/4*sin(circle)*cos(THETHA) + y(n,round(k)); %y-location of particle

            if variables == 3
                color = [abs(v(n,round(k)))/max(max(abs(v))) 0 0];
            else
                color = [abs(v(n,round(k)))/max(max(abs(v))) 0 abs(u(n,round(k)))/max(max(abs(u)))];
            end
            %color is based on the max velocity in v and u direction
            %red = high v, blue = high u, or combo

            p1 = fill(xx,yy,color);
            ellipse_t = fit_ellipse(x_snap, y_snap, fig1);
        end 