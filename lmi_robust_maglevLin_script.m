%% define the important values
f       = 2*pi*1;
A       = [0 1;0 0];
B       = [0; 826];
lambda  = 10;
amp     = 1;

alpha = [0.2314 0.7848 -0.0624 0.0095 0.0215]';
time = [0 1.2];

wn = 60;
zt = 0.707;
Bmul = 1;
Q  = eye(length(A));

%% initialization: calculate K1, K2, K3, Am, Bm and P 
Am = [0 1;-(wn^2) -2*wn*zt];
Bm = [0; wn^2] * Bmul;

% MRAC
P = lyap(Am',Q);

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

numStates = length(A) + length(Am) + length(alpha);

%% solve the system and plot results

[to,xo] = ode15s(@(t,x) nextState(t,x,A,B,K1,K2,K3,Am,Bm,alpha,P,lambda,f), time, zeros(numStates,1));

h  = xo(:,1);
hm = xo(:,3);
s  = xo(:,2);
sm = xo(:,4);
u  = square(f,to);

plot(to, h);
hold on
plot(to, u);
hold off
title('Plant Outputs');
legend('height', 'input');
xlabel('Time (s)');
ylabel('Height (cm)');

figure
plot(to, hm);
hold on
plot(to, u);
hold off
title('Reference Model Outputs');
legend('height', 'input');
xlabel('Time (s)');
ylabel('Height (cm)');

figure
plot(to, h - hm);
hold on
plot(to, s - sm);
title('Errors');
legend('height', 'speed');
xlabel('Time (s)');
ylabel('Height (cm, cm/s)');


%% functions

function dxdt = nextState(t, states,A,B,K1,K2,K3,Am,Bm,alpha,P,lambda,f)
    
    u       = square(f,t);
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
    global amp
    u = amp * (double(sin(f*t) >= 0) - 0.5) * 2;
end

function T = THETA(x)
    T = [x(1); x(2); abs(x(1))*x(2); abs(x(2))*x(2); x(1)^3];
end