global A B K1 K2 K3 Am Bm alpha P lambda


A   = [0 1;0 0];
B   = [0; 1];
% K3  = [-7.49 -6.56];

% the reference model spec
wn = 0.4;
zt = 0.707;

Am = [0 1;-(wn^2) -2*wn*zt];
Bm = [0; wn^2];

% MRAC
Q = eye(length(A));
P = lyap(Am',Q);
lambda = 100;

% calculate K1 and K2 using pole placement
K1 = place(A, B, eig(Am));
K2 = B \ Bm;

% print K1 and K2
display(K1);
display(K2);

% calculate K3 using LMIs
gamma = 1.8;
k = gamma^2;
setlmis([])

% define the LMI variables: X & W

% X is symmetric and same size as Am
X = lmivar(1, [size(Am,1) 1]);

% W is rectangular and size [1 length(Am)]
W = lmivar(2, [1 length(Am)]);

% the first LMI: X > 0
lmiterm([-1 1 1 X],1,1);

% the second LMI: [...] < 0
lmiterm([2 1 1 X], Am, 1, 's');  % Am*X + (Am*X)'
lmiterm([2 1 1 W], B, 1, 's');   % B*W + (B*W)'
lmiterm([2 1 2 0], B);           % B
lmiterm([2 1 3 -X], 1, 1);       % X' = X

lmiterm([2 2 1 0], B');          % B'
lmiterm([2 2 2 0], -1);          % -I
lmiterm([2 2 3 0], 0);           % zero

lmiterm([2 3 1 X], 1, 1);        % X
lmiterm([2 3 2 0], 0);           % 0
lmiterm([2 3 3 0], -k);        % -k*I

LMIs = getlmis;
[TMIN, XFEAS] = feasp(LMIs);
Xs = dec2mat(LMIs, XFEAS, X);
Ws = dec2mat(LMIs, XFEAS, W);
K3 = Ws * inv(Xs);

% the disturbance coeffs
alpha = [0.2314 0.7848 -0.0624 0.0095 0.0215]';

time = [0 120];

numStates = length(A) + length(Am) + length(alpha);

[to,xo] = ode45(@(t,x) nextState(t,x), time, zeros(numStates,1));

plot(to, rad2deg(xo(:,1)));
hold on
plot(to, rad2deg(xo(:,2)));
hold off
title('Plant Outputs');
legend('angle', 'speed');

figure
plot(to, rad2deg(xo(:,3)));
hold on
plot(to, rad2deg(xo(:,4)));
hold off
title('Reference Model Outputs');
legend('angle', 'speed');



function dxdt = nextState(t, states) 
    global A B K1 K2 K3 Am Bm alpha P lambda
    
    u       = square(pi/30,t);
    x       = states(1:2);
    xm      = states(3:4);
    theta   = states(5:end);
    T       = THETA(x);
    G       = lambda * eye(length(T));
    dxdt    = zeros(4 + length(T),1);
    
    % the baseline
    xdot = A*x + B*(K2*u) + B*(alpha' * T) - B*(theta' * T) - B*K1*x;
    xmdot = Am * xm + Bm*u;
    
    % now the adaptive part
    % em = x - xm
    
    % the state error feedback
    em = x - xm;
    xdot = xdot + B*(K3 * em);
    
    % the adaptive law
    % find the rate of change of theta
    thetadot = G * T * em' * P * B;
    
    dxdt(1:2) = xdot;
    dxdt(3:4) = xmdot;
    dxdt(5:end) = thetadot;
end 

function u = square(f, t)
    u = deg2rad(15) * (double(sin(f*t) >= 0) - 0.5) * 2;
end

function T = THETA(x)
    T = [x(1); x(2); abs(x(1))*x(2); abs(x(2))*x(2); x(1)^3];
end

function d = delta(a, x)
    d = sum(a .* THETA(x));
end