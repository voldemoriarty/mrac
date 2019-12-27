m = 0.12; % the mass of the magnet in kg
a = 1.65; % the a coefficient
b = 6.20; % the b coefficient
g = 9.81; % the gravity acceleration

% define the timespan of the simulation
time = [0 10];

% define the step time
step_time = 0;

% define the initial condition
y0 = [0 0];

% define the amplitude of input
% (the control effort)
amp = 3 * 10^4;

% define the ode45 solver
[to,yo] = ode45(@(t,y) next_state(t,y,step_time,m,a,b,g,amp), time, y0);

% plot the results
plot (to, yo(:,1));

function dydt = next_state (t,y,step_time,m,a,b,g,amp)
    dydt = zeros(2,1);
    dydt(1) = y(2);
    u = double(t >= step_time) .* amp;
    dydt(2) = 1/m * ( (u)/(a*((y(1)+b)^4)) - m*g);
end