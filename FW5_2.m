clear all
close all

m = 1; % mass 1 kg
l = 1; % length = 1 m
g = -10; % 10 m/s^2 
M = 10; % cart mass 10 kg
a = 1/(M+m);
dt =  0.1; % Time step for forward euler, seconds
Tend = 20; % Simulation end time, seconds
theta0 = -pi/50; % Initial angle of pendulum

Alin = [0 1 0 0;
     0 0 (-m*a*g)/(2-m*a) 0;
     0 0 0 1;
     0 0 g/(2*l-m*a*l) 0];
Blin = [0; (2*a/(2-m*a)); 0; (-a/(2*l-m*a*l))];
Clin = [1 0 0 0];

syms s
TF = Clin*inv(s*eye(4)-Alin)*Blin;
pretty(simplify(TF));

Pc = [-1 -1+i, -1-i, -2];
Po = [-20 -25 -30 -35];

K = place(Alin, Blin, Pc);
L = place(Alin', Clin', Po);


syms xdot x theta thetadot u;
% introduce state variables:
x1= x; 
x2= xdot; 
x3= theta; 
x4= thetadot;
X = [x1;x2;x3;x4];

A(1,1) = x2;
A(2,1) = (-m*a*g*sin(x3)*cos(x3)+2*a*x4^2*sin(x3)*m*l)/(2-m*a*cos(x3)^2);
A(3,1) = x4;
A(4,1) = (g*sin(x3)-m*a*l*x4^2*sin(x3)*cos(x3)) / (2*l-m*a*l*cos(x3)^2);

B(1,1) = 0*x1;
B(2,1) = 2*a/(2-m*a*cos(x3)^2);
B(3,1) = 0;
B(4,1) = -a*cos(x3)/(2*l-m*a*l*cos(x3)^2);

u = -K*X;
XDOT = A + B*u;
% Initialize values for each state

x = 0;
xdot = 0;
theta = theta0;
thetadot = 0;
out = zeros(Tend/dt,4);
out(1,:) = [x xdot theta thetadot];

for i = 1:1:Tend/dt
    Xd = subs(XDOT);
    x_new = x + Xd(1)*dt;
    xdot_new = xdot + Xd(2)*dt;
    theta_new = theta + Xd(3)*dt;
    thetadot_new = thetadot + Xd(4)*dt;
    x = double(x_new);
    xdot = double(xdot_new);
    theta = double(theta_new);
    thetadot = double(thetadot_new);
    out(round(i),:) = [x xdot theta thetadot];
    %{ 
    %% Add a disturbance
    if i == 300
        thetadot = 0.5;
    end
    %}
end

x_out = out(:,1);
xdot_out = out(:,2);
theta_out = out(:,3);
thetadot_out = out(:,4);

t = dt:dt:Tend;
figure(1)
subplot(2,2,1)
plot(t',x_out)
title('X Position')
xlabel('time (sec)')
ylabel('X distance (m)')
subplot(2,2,2)
plot(t',xdot_out)
title('X Velocity')
xlabel('time (sec)')
ylabel('X Velocity (m/s)')
subplot(2,2,3)
plot(t',theta_out)
title('Theta Angle')
xlabel('time (sec)')
ylabel('Theta (rad)')
subplot(2,2,4)
plot(t',thetadot_out)
title('Thetadot Angular Velocity')
xlabel('time (sec)')
ylabel('Thetadot (rad/sec)')

xdim = ceil(max(max(x_out),abs(min(x_out)))+0.5);
figure(2)
xlim([-xdim,xdim])
ylim([-0.1,1.5])
daspect([1 1 1])
ground = line([-xdim xdim],[0 0],'color','k');
video = VideoWriter('ECE602_FW5_part1');
set(video,'FrameRate',(30));
open(video);

for i = 1:1:round(Tend/dt)
    mpos = [l*sin(theta_out(i)), l*cos(theta_out(i))];
    cart = rectangle('Position',[x_out(i)-0.25, 0, 0.5, 0.25], 'FaceColor', 'b');
    cartx = get(cart, 'Position');
    bar = line('xdata',[cartx(1)+0.25, cartx(1)+0.25+mpos(1)],'ydata', [0.25,0.25+mpos(2)],'Linewidth',3);
    drawnow
    Frame = getframe;
    writeVideo(video,Frame);
    delete(cart)
    delete(bar)
    
end
close(video)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms y yhat
x = 1;
xdot = 0;
theta = theta0;
thetadot = 0;
out = zeros(Tend/dt,4);
out(1,:) = [x xdot theta thetadot];

Xhat = X;
y =  Clin*X;
yhat = Clin*Xhat;

u = -K*Xhat;
%XDOT = A + B*u;
%Xdh = 

for i = 1:1:Tend/dt
    Xd = subs(A + B*u);
    Xdhat = subs(A + B*u+L*(y-yhat));
    
    x_new = x + Xd(1)*dt;
    xdot_new = xdot + Xd(2)*dt;
    theta_new = theta + Xd(3)*dt;
    thetadot_new = thetadot + Xd(4)*dt;
    x = double(x_new);
    xdot = double(xdot_new);
    theta = double(theta_new);
    thetadot = double(thetadot_new);
    
    Xhat_new = double(subs(Xhat)) + double(subs(Xdhat))*dt;
    Xhat = double(subs(Xhat_new));
    u = -K*Xhat;
    y_new = Clin*double(subs(X));
    yhat_new = Clin*double(subs(Xhat));
    out(round(i),:) = [x xdot theta thetadot];
end

x_out = out(:,1);
xdot_out = out(:,2);
theta_out = out(:,3);
thetadot_out = out(:,4);

%% Plot
t = dt:dt:Tend;
figure(3)
subplot(2,2,1)
plot(t',x_out)
title('X Position')
xlabel('time (sec)')
ylabel('X distance (m)')
subplot(2,2,2)
plot(t',xdot_out)
title('X Velocity')
xlabel('time (sec)')
ylabel('X Velocity (m/s)')
subplot(2,2,3)
plot(t',theta_out)
title('Theta Angle')
xlabel('time (sec)')
ylabel('Theta (rad)')
subplot(2,2,4)
plot(t',thetadot_out)
title('Thetadot Angular Velocity')
xlabel('time (sec)')
ylabel('Thetadot (rad/sec)')

%% Simulate
xdim = ceil(max(max(x_out),abs(min(x_out)))+0.5);
figure(4)
xlim([-xdim,xdim])
ylim([-0.1,1.5])
daspect([1 1 1])
ground = line([-xdim xdim],[0 0],'color','k');
video = VideoWriter('ECE602_FW5_part2');
set(video,'FrameRate',(30));
open(video);

for i = 1:1:round(Tend/dt)
    mpos = [l*sin(theta_out(i)), l*cos(theta_out(i))];
    cart = rectangle('Position',[x_out(i), 0, 0.5, 0.25], 'FaceColor', 'b');
    cartx = get(cart, 'Position');
    bar = line('xdata',[cartx(1) + 0.25, cartx(1)+0.25+mpos(1)],'ydata', [0.25,0.25+mpos(2)],'Linewidth',3);
    drawnow
    Frame = getframe;
    writeVideo(video,Frame);
    delete(cart)
    delete(bar)
    
end
close(video)
%}