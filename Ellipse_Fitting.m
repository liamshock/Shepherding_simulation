% Ellipse from data

% check a single snapshot 
tmidx = 50;
%x_curr = x(:,tmidx); y_curr = y(:,tmidx);
x_curr = x(:,2); y_curr = y(:,2);
pos_curr = [x_curr';y_curr'];
% Fit ellipse to data 
[z, a, b, alpha] = fitellipse(pos_curr);
Q = [cos(alpha), -sin(alpha); 
     sin(alpha), cos(alpha)];
theta = linspace(0,2*pi,360);
PosEllip = z + Q*[a * cos(theta); b * sin(theta)];
Orient = mod(alpha*180/pi,360);
Eccen = b\a;

K = convhull(x(:,tmidx),y(:,tmidx));
figure 
plot(x_curr,y_curr,'.k'); axis equal tight 
hold on 
plot(x_curr(K),y_curr(K),'b','linewidth',3);
hold on
plot(PosEllip(1,:),PosEllip(2,:),'r','linewidth',3);
