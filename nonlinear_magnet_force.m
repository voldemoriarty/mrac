m = 0.12; % the mass of the magnet in kg
a = 1.65; % the a coefficient
b = 6.20; % the b coefficient
g = 9.81; % the gravity acceleration in m/s/s

n = 1000;

% plot the characteristic for control effort against position
position = linspace(0,5,n);

coeffort = (m*g*a) .* ((position + b) .^ 4);

plot(position, coeffort);

title('u_1 actuator characteristic');
xlabel('Magnet Position (cm)');
ylabel('Control Effort (u_1)');