% Enhanced 3D Rocket Trajectory Simulation

% Constants
g = 9.81;              % Acceleration due to gravity (m/s^2)
m0 = 100;              % Initial mass of the rocket (kg)
fuel_rate = 1;         % Fuel consumption rate (kg/s)
Thrust_mag = 1000;     % Magnitude of thrust (N)
Cd = 0.2;              % Drag coefficient
A = 0.1;               % Cross-sectional area of the rocket (m^2)
rho = 1.2;             % Air density (kg/m^3)
wind = [1, 1, 0];      % Wind velocity (m/s)
thrust_direction = [0, 0, 1]; % Thrust direction (unit vector)

% Initial conditions
position0 = [0, 0, 0];  % Initial position (m)
velocity0 = [0, 0, 0];  % Initial velocity (m/s)
dt = 0.1;              % Time step (s)
T = 100;                % Total simulation time (s)
time = 0:dt:T;         % Time vector
N = length(time);

% Variables to hold the simulation results
position = zeros(N, 3);
velocity = zeros(N, 3);
mass = zeros(N, 1);
position(1,:) = position0;
velocity(1,:) = velocity0;
mass(1) = m0;

% Main simulation loop
for i = 2:N
    m = mass(i-1) - fuel_rate * dt; % Updated mass
    if m < m0/2 % Stop thrust after burning half the fuel
        Thrust = [0, 0, 0];
    else
        Thrust = Thrust_mag * thrust_direction; % Thrust force
    end

    % Relative velocity (rocket's velocity minus wind)
    relative_velocity = velocity(i-1,:) - wind;
    
    % Drag force
    Drag = 0.5 * Cd * A * rho * norm(relative_velocity) * relative_velocity;
    
    % Net force (Thrust minus drag and gravity)
    F_net = Thrust - Drag - [0, 0, m * g];
    
    % Acceleration
    a = F_net / m;
    
    % Update velocity and position using simple Euler integration
    velocity(i,:) = velocity(i-1,:) + a * dt;
    position(i,:) = position(i-1,:) + velocity(i-1,:) * dt;
    
    mass(i) = m; % Store mass
    
    % Stop simulation if rocket hits the ground
    if position(i,3) < 0
        position(i,3) = 0;
        velocity(i,:) = [0, 0, 0];
        break;
    end
end

% Plotting the results
figure
subplot(3,1,1);
plot3(position(:,1), position(:,2), position(:,3))
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')
title('3D Rocket Trajectory')
grid on

subplot(3,1,2);
plot(time, position)
xlabel('Time (s)')
ylabel('Position (m)')
legend('X', 'Y', 'Z')
title('Rocket Position vs Time')
grid on

subplot(3,1,3);
plot(time, mass)
xlabel('Time (s)')
ylabel('Mass (kg)')
title('Rocket Mass vs Time')
grid on
