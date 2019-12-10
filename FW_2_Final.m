%FW_2 Final
clc
clear all
close all

L1 = 0.5; % Length of the first pendulum (m)
L2 = 0.5; % Length of the second pendulum (m)
m1 = 1; % Mass 1 (kg)
m2 = 1; % Mass 2 (kg)
g = 9.8; % acceleration of gravity (m/s^2)

% Initial conditions for theta and theta dot
theta_init = [135 0 180 0]; % [theta2, theta2_dot, theta1, theta1_dot]
theta_init_rad = theta_init.*pi()./180; % Covert inital position to radians
T_end = 15*60;
t_i = 0.0025;
i_end = T_end / t_i + 1;
t_int = 0:t_i:T_end;

[t,y] = ode45(@pend, t_int, theta_init_rad);

theta_1 = y(:,1)';
theta_2 = y(:,3)';

% Using the original position calculation with the theta vectors
m1_posx = L1*sin(theta_1);
m1_posy = -L1*cos(theta_1);
m2_posx = L1*sin(theta_1) + L2*sin(theta_2);
m2_posy = -L1*cos(theta_1) - L2*cos(theta_2);

figure(1)
hold on
axis([-1.1 1.1 -1.1 1.1])
p = zeros(1,4);
p(1) = plot(m1_posx, m1_posy,'k');
p(2) = plot(m2_posx, m2_posy,'m');

origin = [0, 0];
m1 = [L1*sind(theta_init(1)), -L1*cosd(theta_init(1))];
m2 = [m1(1) + L2*sind(theta_init(3)), m1(2) -  L2*cosd(theta_init(3))];
p(3) = line([origin(1),m1(1)], [origin(2), m1(2)],'Color','k','LineWidth',4);
p(4) = line([m1(1),m2(1)], [m1(2), m2(2)],'Color','b','LineWidth',4);

legend([p(1) p(2)],'Mass 1 Position', 'Mass 2 Position')
title('Double Pendulum Simulation')
xlabel('x position (m)')
ylabel('y position (m)')
hold off

figure(2)
hold on
axis([-1 1 -1 1])
plot(0,0,'.','markersize',20,'color','k')
title('Double Pendulum Simulation')
xlabel('x position (m)')
ylabel('y position (m)')
q = zeros(1,8);
F(1,i_end-1) = struct('cdata','','colormap','');
for i = 1:i_end-1
    m1 = [L1*sin(theta_1(i)), -L1*cos(theta_1(i))];
    m2 = [m1(1) + L2*sin(theta_2(i)), m1(2) -  L2*cos(theta_2(i))];
    if mod(i,2) == 1
        q(1) = line([origin(1),m1(1)], [origin(2), m1(2)],'Color','k','LineWidth',4)
        q(2) = line([m1(1),m2(1)], [m1(2), m2(2)],'Color','b','LineWidth',4)
        q(3) = plot(m1(1),m1(2),'.','markersize',20,'color', 'm')
        q(4) = plot(m2(1),m2(2),'.','markersize',20,'color', 'c')
        if i > 1
            delete(q(5))
        	delete(q(6))
            delete(q(7))
        	delete(q(8))
        end
    end
    if mod(i,2) == 0
        q(5) = line([origin(1),m1(1)], [origin(2), m1(2)],'Color','k','LineWidth',4)
        q(6) = line([m1(1),m2(1)], [m1(2), m2(2)],'Color','b','LineWidth',4)
        q(7) = plot(m1(1),m1(2),'.','markersize',20,'color', 'm')
        q(8) = plot(m2(1),m2(2),'.','markersize',20,'color', 'c')
        delete(q(1))
        delete(q(2))
        delete(q(3))
        delete(q(4))
    end

    %pause(t_i)
    F(i) = getframe(gcf);
end

video = VideoWriter('mezmerize.avi');
video.FrameRate = 1/t_i;
open(video);
for j=1:length(F)
    frame = F(j);
    writeVideo(video,frame);
end
close(video);

function dx = pend(t, x)
    L1 = 0.5; % Length of the first pendulum (m)
    L2 = 0.5; % Length of the second pendulum (m)
    m1 = 1; % Mass 1 (kg)
    m2 = 1; % Mass 2 (kg)
    t = 0.1;
    g = 9.8; % acceleration of gravity (m/s^2)
    
    dx = zeros(4, 1);
    dx(1) = x(2);
    dx(2) = (m2*cos(x(1) - x(3))*(- L1*sin(x(1) - x(3))*x(2)^2 + g*sin(x(3))))/(- L1*m2*cos(x(1) - x(3))^2 + L1*m1 + L1*m2) - (L2*m2*sin(x(1) - x(3))*x(4)^2 + g*sin(x(1))*(m1 + m2))/(- L1*m2*cos(x(1) - x(3))^2 + L1*m1 + L1*m2);
    dx(3) = x(4);
    dx(4) = (cos(x(1) - x(3))*(L2*m2*sin(x(1) - x(3))*x(4)^2 + g*sin(x(1))*(m1 + m2)))/(- L2*m2*cos(x(1) - x(3))^2 + L2*m1 + L2*m2) - ((m1 + m2)*(- L1*sin(x(1) - x(3))*x(2)^2 + g*sin(x(3))))/(- L2*m2*cos(x(1) - x(3))^2 + L2*m1 + L2*m2);

end