clear all
clc
close all
addpath('C:/Users/sabin/Desktop/KoopmanMPC-master/KoopmanMPC-master/Resources')
addpath('C:/Users/sabin/Desktop/KoopmanMPC-master/KoopmanMPC-master/Resources/qpOASES-3.1.0/interfaces/matlab') 

rng(2141444)



%% *************************** Dynamics ***********************************

f_u = @pend_dynamics;
n = 4;
m = 1; % number of control inputs

xrange = 1;
urange = 20;
%% ************************** Discretization ******************************

deltaT = 0.05;
%Runge-Kutta 4
k1 = @(t,x,u) (  f_u(t,x,u) );
k2 = @(t,x,u) ( f_u(t,x + k1(t,x,u)*deltaT/2,u) );
k3 = @(t,x,u) ( f_u(t,x + k2(t,x,u)*deltaT/2,u) );
k4 = @(t,x,u) ( f_u(t,x + k1(t,x,u)*deltaT,u) );
f_ud = @(t,x,u) ( x + (deltaT/6) * ( k1(t,x,u) + 2*k2(t,x,u) + 2*k3(t,x,u) + k4(t,x,u)  )   );

%% ************************** Collect data ********************************
tic
disp('Starting data collection')
Nsim = 200;
Ntraj = 1000;
Cy = [1 0 0 0; 0 0 1 0]; % Output matrix: y = Cy*x
ny = size(Cy,1);

% Random forcing
Ubig = 2*urange*rand([Nsim Ntraj]) - urange;

% Random initial conditions
Xcurrent = (rand(n,Ntraj)*xrange*2 - xrange);
MK = Xcurrent;
X = []; Y = []; U = [];
%%
for i = 1:Nsim
    Xnext = f_ud(0,Xcurrent,Ubig(i,:));
    X = [X Xcurrent];
    Y = [Y Xnext];
    U = [U Ubig(i,:)];
    Xcurrent = Xnext;
end
fprintf('Data collection DONE, time = %1.2f s \n', toc);


%% ************************** Basis functions *****************************
n_zeta = ny;
basisFunction = 'rbf';
Nrbf = 100;
cent = rand(n,Nrbf)*xrange*2 - xrange; % RBF centers
rbf_type = 'thinplate';
theta_max = pi;
liftFun = @(xx)( [xx;rbf(xx,cent,rbf_type)] );
Nlift = Nrbf + n;

%%
figure
plot([0:Nsim-1]*deltaT, Y(1,1+10000:Nsim+10000),'-r','linewidth',2);

%% ******************************* Lift ***********************************
disp('Starting LIFTING')
tic
Xlift = liftFun(X);
Ylift = liftFun(Y);
fprintf('Lifting DONE, time = %1.2f s \n', toc);

%% ********************** Build predictor *********************************

disp('Starting REGRESSION for A,B,C')
tic
W = [Ylift ; X];
V = [Xlift ; U];
VVt = V*V';
WVt = W*V';
ABC = WVt * pinv(VVt);
Alift = ABC(1:Nlift,1:Nlift);
Blift = ABC(1:Nlift,Nlift+1:end);
Clift = ABC(Nlift+1:end,1:Nlift);
fprintf('Regression for A, B, C DONE \n');

% Residual
fprintf( 'Regression residual %f \n', norm(Ylift - Alift*Xlift - Blift*U,'fro') / norm(Ylift,'fro') );

%% ************************* Predictor comparison *************************
Tmax = 5;
Nsim = Tmax/deltaT;
uprbs = (2*urange*myprbs(Nsim,0.5) - urange);
u_dt = @(i)(  uprbs(i+1) );
%u_dt = @(i)((-1).^(round(i/30))); % control signal

% Initial condition
x0 = rand(n,1)-0.5;

%x0 = [ 5.0; 0.; pi/2-.15; 0. ];
x_true = x0;

% Lifted initial condition
xlift = liftFun(x0);

% Local linearization predictor at x0
x = sym('x',[n;1]); u = sym('u',[1;1]);
Ac_x0 = double(subs(jacobian(f_u(0,x,u),x),[x;u],[x0;0]));
Bc_x0 = double(subs(jacobian(f_u(0,x,u),u),[x;u],[x0;0]));
c_x0 = double(subs(f_u(0,x,u),[x;u],[x0;0])) - Ac_x0*x0 - Bc_x0*0;
ABc = expm([Ac_x0 [Bc_x0 c_x0] ; zeros(2,6)]*deltaT); % discretize
Ad_x0 = ABc(1:n,1:n); Bd_x0 = ABc(1:n,3); cd_x0 = ABc(1:n,n);
X_loc_x0 = x0;

% Local linearization predictor at 0
x = sym('x',[n;1]); u = sym('u',[1;1]);
Ac_0 = double(subs(jacobian(f_u(0,x,u),x),[x;u],[0;0;0;0;0]));
Bc_0 = double(subs(jacobian(f_u(0,x,u),u),[x;u],[0;0;0;0;0]));
c_0 = double(subs(f_u(0,x,u),[x;u],[0;0;0;0;0])) - Ac_0*[0;0;0;0] - Bc_0*0;
ABc = expm([Ac_0 [Bc_0 c_0] ; zeros(2,6)]*deltaT); % discretize
Ad_0 = ABc(1:n,1:n); Bd_0 = ABc(1:n,3); cd_0 = ABc(1:n,n); 
X_loc_0 = x0;


% Simulate
for i = 0:Nsim-1
    % Koopman predictor
    xlift = [xlift, Alift*xlift(:,end) + Blift*u_dt(i)]; % Lifted dynamics
    
    % True dynamics
    x_true = [x_true, f_ud(0,x_true(:,end),u_dt(i)) ];
    
    % Local linearization predictor at x0
 

    X_loc_x0 = [X_loc_x0, Ad_x0*X_loc_x0(:,end) + Bd_x0*u_dt(i) + cd_x0];
    
    % Local linearization predictor at 0
    X_loc_0 = [X_loc_0, Ad_0*X_loc_0(:,end) + Bd_0*u_dt(i) + c_0];
    
end

x_koop = Clift * xlift; % Koopman predictions

%%

lw = 4;

figure
plot([0:Nsim-1]*deltaT,u_dt(0:Nsim-1),'linewidth',lw); hold on
title('Control input $u$', 'interpreter','latex'); xlabel('Time [s]','interpreter','latex');
set(gca,'fontsize',20)

figure
plot([0:Nsim]*deltaT,x_true(3,:),'linewidth',lw); hold on
plot([0:Nsim]*deltaT,x_koop(3,:), '--r','linewidth',lw)
plot([0:Nsim]*deltaT,X_loc_x0(3,:), '--g','linewidth',lw-1)
plot([0:Nsim]*deltaT,X_loc_0(3,:), '--k','linewidth',lw-1)
axis([0 Tmax min(x_koop(3,:))-0.15 max(x_koop(3,:))+0.15])
title('Predictor comparison - $x_2$','interpreter','latex'); xlabel('Time [s]','interpreter','latex');
set(gca,'fontsize',20)
LEG = legend('True','Koopman','Local at $x_0$','Local at 0','location','southwest');
set(LEG,'interpreter','latex')

figure
plot([0:Nsim]*deltaT,x_true(1,:),'linewidth',lw); hold on
plot([0:Nsim]*deltaT,x_koop(1,:), '--r','linewidth',lw)
plot([0:Nsim]*deltaT,X_loc_x0(1,:), '--g','linewidth',lw-1)
plot([0:Nsim]*deltaT,X_loc_0(1,:), '--k','linewidth',lw-1)
axis([0 Tmax min(x_koop(1,:))-0.1 max(x_koop(1,:))+0.1])
title('Predictor comparison - $x_1$','interpreter','latex'); xlabel('Time [s]','interpreter','latex');
set(gca,'fontsize',20)
LEG = legend('True','Koopman','Local at $x_0$','Local at 0','location','southwest');
set(LEG,'interpreter','latex')



%% ********************** Feedback control ********************************
disp('Press any key for feedback control')
pause

Tmax = 5; % Simlation legth
Nsim = Tmax/deltaT;

ymin = -5;
ymax = 5;

umin = -10;
umax = 10;

thetamin = -2*pi;
thetamax = 2*pi;

%x0 = [ 0.0; 0.; pi/2-.15; 0. ];
x0 = rand(n,1)-0.5;
%x0 = [ 5.0; 0.; pi/2-.15; 0. ];
% Define Koopman controller
C = zeros(2,Nlift); C(1,1) = 1; C(2,3) = 1;

% Weight matrices
Q = 100;
R = 1;
% Prediction horizon
Tpred = 1;
Np = round(Tpred / deltaT);
% Constraints
xlift_min = [ymin ; nan; -10; nan(Nlift-3,1)];
xlift_max = [ymax ; nan; 10; nan(Nlift-3,1)];

% Build Koopman MPC controller
koopmanMPC  = getMPC(Alift,Blift,C,0,Q,R,Q,Np,umin, umax, xlift_min, xlift_max,'qpoases');

% Initial condition for the delay-embedded state (assuming zero control in the past)
x = x0;
zeta0 = x;


x_koop = x0; x_loc = x0;
zeta = zeta0; 

XX_koop = x0; UU_koop = [];
XX_loc = x0; UU_loc = [];

% Get Jacobian of the true dynamics (for local linearization MPC)
x = sym('x',[n 1]); syms u;
f_ud_sym = f_ud(0,x,u);
u_loc = 0;

Jx = jacobian(f_ud_sym,x);
Ju = jacobian(f_ud_sym,u);

wasinfeas= 0;
ind_inf = [];

% Closed-loop simultion start
UU_koop1 = [];

yrr1 = X(1,1:Nsim);
yrr2 = X(3,1:Nsim);

yrr = 0.5*cos(2*pi*[1:Nsim] / Nsim);
yrrc = [yrr1;yrr];
for i = 0:Nsim-1
    if(mod(i,10) == 0)
        fprintf('Closed-loop simulation: iterate %i out of %i \n', i, Nsim)
    end
    
    % Current value of the reference signal
    yr = yrrc(:,i+1);

    % Koopman MPC
    xlift = liftFun(zeta); % Lift
    u_koop = koopmanMPC(xlift, [0;pi]); % Get control input
    x_koop = f_ud(0,x_koop,u_koop); % Update true state
    zeta = x_koop;  % Update delay-embedded state
    
    % Local linearization MPC
    Aloc = double(subs(Jx,[x;u],[x_loc;u_loc])); % Get local linearization
    Bloc = double(subs(Ju,[x;u],[x_loc;u_loc]));
    cloc = double(subs(f_ud_sym,[x;u],[x_loc;u_loc])) - Aloc*x_loc - Bloc*u_loc;
    [U_loc,~,optval] = solveMPCprob(Aloc,Bloc,Cy,cloc,Q,R,Q,Np,umin, umax,[ymin;nan;-10;nan],[ymax;nan;10;nan],x_loc,yr); % Get control input
    u_loc = U_loc(1:m,1);
    if(optval == Inf) % Detect infeasibility
        ind_inf = [ind_inf i];
        wasinfeas = 1;
    end
    x_loc = f_ud(0,x_loc,u_loc); % Update true state
    
    % Store values
    XX_koop = [XX_koop x_koop];
    UU_koop = [UU_koop u_koop];
    XX_loc = [XX_loc x_loc];
    UU_loc = [UU_loc u_loc];
end

if(isempty(ind_inf))
    ind_inf = Nsim;
end

%% Plot (feedback control)
lw_koop = 3;


ww = repmat(yr, Nsim,1);
ww = ww(1:Nsim, :);
% Control signal
figure
p3 = plot([0:Nsim]*deltaT,umax*ones(Nsim+1,1),'-k','linewidth',lw_koop-1); hold on
p4 = plot([0:Nsim]*deltaT,umin*ones(Nsim+1,1),'-k','linewidth',lw_koop-1); hold on
p1 = plot([0:ind_inf(1)-1]*deltaT,UU_loc(1:ind_inf(1)),'--g','linewidth',lw_koop); hold on
p2 = plot([0:Nsim-1]*deltaT,UU_koop,'-b','linewidth',lw_koop); hold on
axis([0,Tmax,min(UU_loc) - abs(min(UU_loc))*0.1, max(UU_loc)*1.1] )
LEG  = legend([p2,p1,p3],'K-MPC','L-MPC','Constraint');
set(LEG,'Interpreter','latex','location','southeast')
set(gca,'FontSize',31);
 
figure
% Output (y = x_2)
p1=plot([0:Nsim]*deltaT,ymax*ones(Nsim+1,1),'-k','linewidth',lw_koop-1); hold on
p2=plot([0:Nsim]*deltaT,ymin*ones(Nsim+1,1),'-k','linewidth',lw_koop-1);
p3=plot([0:Nsim]*deltaT,XX_koop(1,:),'-b','linewidth',lw_koop); hold on
p4=plot([0:ind_inf(1)-1]*deltaT,XX_loc(1,1:ind_inf(1)),'--g','linewidth',lw_koop);
p5=plot([0:Nsim-1]*deltaT, yrr1,'--r','linewidth',lw_koop);
LEG  = legend([p3,p4,p5,p2],'K-MPC','L-MPC','Reference','Constraints');
set(LEG,'Interpreter','latex','location','southeast')
set(LEG,'Fontsize',18)
axis([0,Tmax,-5,5])
set(gca,'FontSize',20);

%%
% 
% figure
% plot([0:Nsim-1]*deltaT, X(3,1+500:Nsim+500),'-r','linewidth',lw_koop);
% %%
% yrrx = 2*cos(2*pi*[1:Nsim] / Nsim);
% figure
% plot([0:Nsim-1]*deltaT, yrrx,'-r','linewidth',lw_koop);
%%
function dYdt = pend_dynamics(t,y,u)
    g = 9.8; 
    L = 1.5; 

    m = 1.0; %mass of bob (kg)
    M = 5.0;  % mass of cart (kg)
    
    d1 = 1.0;
    d2 = 0.5;

  
    x_ddot = u - (m*L).*y(4,:).*y(4,:).*cos(y(3,:)) + (m*g).*cos(y(3,:)).*sin(y(3,:));

    x_ddot = x_ddot ./ ( (M+m)-m.* sin(y(3,:)).*sin(y(3,:)) );

    theta_ddot = -g/L .* cos(y(3,:)) - 1./L .* sin( y(3,:) ) .* x_ddot;
    
    damping_x =  - d1.*y(2,:);
    damping_theta =  - d2.*y(4,:);
    
    dYdt = [ y(2,:); x_ddot+damping_x; y(4,:); theta_ddot+damping_theta];
    
end




