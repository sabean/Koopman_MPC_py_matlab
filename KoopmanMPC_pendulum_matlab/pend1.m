clear all
clc
close all

addpath('C:/Users/sabin/Desktop/KoopmanMPC-master/KoopmanMPC-master/Resources')
addpath('C:/Users/sabin/Desktop/KoopmanMPC-master/KoopmanMPC-master/Resources/qpOASES-3.1.0/interfaces/matlab') 


%% ****************************** Dynamics ********************************
%rng(115123)

n = 4; % Number of states
m = 1; % Number of control inputs
f_u = @pend; % Dynamics

umin = -1;
umax = 1;


x1min = -2.5;
x1max = 0.15;

xdotmin = -0.03;
xdotmax = 0.03;

thetamin = -0.89;
thetamax = 1.9;

thetadotmin = -1;
thetadotmax = 1;

uconstraints =[umin;umax];

stateconstrainsmin = [x1min;xdotmin;thetamin;thetadotmin];
stateconstrainsmax = [x1max;xdotmax;thetamax;thetadotmax];



% Discretize
deltaT = 0.01;
%Runge-Kutta 4
k1 = @(t,x,u) (  f_u(t,x,u) );
k2 = @(t,x,u) ( f_u(t,x + k1(t,x,u)*deltaT/2,u) );
k3 = @(t,x,u) ( f_u(t,x + k2(t,x,u)*deltaT/2,u) );
k4 = @(t,x,u) ( f_u(t,x + k1(t,x,u)*deltaT,u) );
f_ud = @(t,x,u) ( x + (deltaT/6) * ( k1(t,x,u) + 2*k2(t,x,u) + 2*k3(t,x,u) + k4(t,x,u)  )   );


%% Collect data
rng(115123)
disp('Starting data collection')
Ntraj = 20; Nsim = 1000;

Cy = [1 0 0 0; 0 0 1 0]; % Output matrix: y = Cy*x
nD = 1; % Number of delays
ny = size(Cy,1); % Number of outputs
n_zeta = (nD+1)*ny + nD*m; % dimension of the delay-embedded "state"

[X,Y,U] = CollectDataPendDelay(deltaT,uconstraints,stateconstrainsmin, stateconstrainsmax,Ntraj,Nsim,n);
fprintf('Data collection DONE \n');

%% Basis functions
basisFunction = 'rbf';
Nrbf = 100;
cent = rand(n_zeta,Nrbf)*2 - 1; % RBF centers
rbf_type = 'polyharmonic';
theta_max = pi;
liftFun = @(xx)( [xx;rbf(xx,cent,rbf_type)] );
Nlift = Nrbf + n_zeta;
%%
figure
NS=50;
plot([0:NS-1]*deltaT, X(1,1+100:NS+100),'-r','linewidth',2);

%% Lift
disp('Starting LIFTING')

Xlift = liftFun(X);
Ylift = liftFun(Y);

%% Regression

disp('Starting REGRESSION for A,B,C')

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
Tmax = 10;
Nsim = Tmax/deltaT;

uprbs = (2*myprbs(Nsim,0.5) - 1);
u_dt = @(i)(  uprbs(i+1) );
%u_dt = @(i)((-1).^(round(i/30)));
f_cont_d = @(t,xx)( f_ud(t,xx,u_dt(t) ));

x0 = rand(n,1)-0.5;
x = x0;

% Delayed initial condition (assume random control input in the past)
xstart = [Cy*x ; NaN(nD*(ny+m),1)];
for i = 1:nD
    urand = 2*rand(m,1) - 1;
    xp = f_ud(0,x,urand);
    xstart = [Cy*xp ; urand; xstart(1:end-ny-m)];
    x = xp;
end

% Local linearization
% xloc = xp;
% x = sym('x',[n 1]); syms u;
% Aloc = double(subs(jacobian(f_ud(0,x,u,deltaT,size(x)),x),[x;u],[xloc;urand]));
% Bloc = double(subs(jacobian(f_ud(0,x,u, deltaT,size(x)),u),[x;u],[xloc;urand]));
% cloc = double(subs(f_ud(0,x,u,deltaT,size(x)),[x;u],[xloc;urand])) - Aloc*xloc - Bloc*urand;

% Inital conditions
x_true = xp;
xlift = liftFun(xstart);

% Simulation
for i = 0:Nsim-1
    
    % True dynamics
    x_true = [x_true, f_ud(0,x_true(:,end),u_dt(i)) ];
    
    % Koopman predictor
    xlift = [xlift Alift*xlift(:,end) + Blift*u_dt(i)];
%     
%     % Local linearization predictor
%     xloc = [xloc Aloc*xloc(:,end) + Bloc*u_dt(i) + cloc];
end


figure
stairs([0:Nsim-1]*deltaT,u_dt(0:Nsim-1),'linewidth',2); hold on
title('Control input'); xlabel('time [s]')


lw_koop = 3;
s11= Cy*x_true;
sc11= Clift(1,:)*xlift;
sc12 = Clift(2,:)*xlift;
% s13= Cy*xloc;

figure

plot([0:Nsim]*deltaT,s11(1,:),'--b','linewidth', lw_koop); hold on
plot([0:Nsim]*deltaT,sc11, '--r','linewidth',lw_koop)
% plot([0:Nsim]*deltaT,s13(1,:), '--g','linewidth',lw_koop-1)

LEG = legend('True','Koopman','Local at $x_0$');
set(LEG,'Interpreter','latex','location','northeast','fontsize',30)
set(gca,'FontSize',25);
axis([0 1 -1.3 0.5])


figure

plot([0:Nsim]*deltaT,s11(2,:),'--b','linewidth', lw_koop); hold on
plot([0:Nsim]*deltaT,sc12, '--r','linewidth',lw_koop)
% plot([0:Nsim]*deltaT,s13(2,:), '--g','linewidth',lw_koop-1)
%plot([0:Nsim]*deltaT,Cy*x_true,'--b','linewidth', lw_koop); hold on
%plot([0:Nsim]*deltaT,Clift(1,:)*xlift, '--r','linewidth',lw_koop)
%plot([0:Nsim]*deltaT,Cy*xloc, '--g','linewidth',lw_koop-1)

LEG = legend('True','Koopman','Local at $x_0$');
set(LEG,'Interpreter','latex','location','northeast','fontsize',30)
set(gca,'FontSize',25);
axis([0 1 -1.3 0.5])



%% ********************** Feedback control ********************************
disp('Press any key for feedback control')
pause

Tmax = 10; % Simlation legth
Nsim = Tmax/deltaT;

ymin = -10;
ymax = 10;
x0 = [ 0.0; 0.; pi; 0. ];
%yrr = 0.5*cos(2*pi*[1:Nsim] / Nsim); % reference


% Define Koopman controller
C = zeros(2,Nlift); C(1,1) = 1; C(2,3) = 1;
% Weight matrices
Q = 200;
R = 1;
% Prediction horizon
Tpred = 1;
Np = round(Tpred / deltaT);
% Constraints
xlift_min = nan(Nlift,1);
xlift_max = nan(Nlift,1);

% Build Koopman MPC controller
koopmanMPC  = getMPC(Alift,Blift,C,0,Q,R,Q,Np,-1, 1, xlift_min, xlift_max,'qpoases');

REF = 'step'; % 'step' or 'cos'
switch REF
    case 'step'
        ymin = -0.6;
        ymax = 0.6;
        yrr11 = 0.3*( -1 + 2*([1:Nsim] > Nsim/2)  ); % reference
        yrr12 = 0.3*( -1 + 2*([1:Nsim] > Nsim/2)  ); % reference
        
    case 'cos'
        ymin = 2.6;ymax = 3.6;
        ymin1 = 0;ymax1 =-8.3 ;
        yrr11 = cos(2*pi*[1:Nsim] / Nsim)*(ymax1-ymin1)+ymin1; % reference
        yrr12 = cos(2*pi*[1:Nsim] / Nsim)*(ymax-ymin)+ymin; % reference
    case 'reference'
        yrr11 = X(1,1:Nsim);
        yrr12 = X(3,1:Nsim); 
end
yrrc = [yrr11;yrr11;yrr12;yrr12];
x0=yrrc(:,1);


% Initial condition for the delay-embedded state (assuming zero control in the past)
x = x0;
zeta0 = [Cy*x ; NaN(nD*(ny+m),1)];
for i = 1:nD
    upast = zeros(m,1);
    xp = f_ud(0,x,upast);
    zeta0 = [Cy*xp ; upast ; zeta0(1:end-ny-m)]; 
    x = xp;
end
x0 = x;

x_koop = x0; x_loc = x0;
zeta = zeta0; % Delay-embedded "state"

XX_koop = x0; UU_koop = [];
% XX_loc = x0; UU_loc = [];

% Get Jacobian of the true dynamics (for local linearization MPC)
% x = sym('x',[n 1]); syms u;
% f_ud_sym = f_ud(0,x,u,deltaT, size(x));
% u_loc = 0;
% Jx = jacobian(f_ud_sym,x);
% Ju = jacobian(f_ud_sym,u);

wasinfeas= 0;
ind_inf = [];

% Closed-loop simultion start
UU_koop1 = [];

yrr1 = X(1,1:Nsim);
yrr2 = X(3,1:Nsim);
% yrrc = [yrr1;yrr2];
yrrc = X(:, 1:Nsim);
yrliftstack = [];
for i = 0:Nsim-1
    if(mod(i,10) == 0)
        fprintf('Closed-loop simulation: iterate %i out of %i \n', i, Nsim)
    end
    
    % Current value of the reference signal
    %yr = yrr(i+1);
    
    yr = yrrc(:,i+1);
    
    % Koopman MPC
    xlift = liftFun(zeta); % Lift
    
    
    yrrlift = [yr(1); yr(3)];
    yrliftstack = [yrliftstack yrrlift];
    u_koop = koopmanMPC(xlift,yrrlift); % Get control input
    x_koop = f_ud(0,x_koop,u_koop); % Update true state
    zeta = [ Cy*x_koop ; u_koop; zeta(1:end-ny-m)]; % Update delay-embedded state
    
%     % Local linearization MPC
%     Aloc = double(subs(Jx,[x;u],[x_loc;u_loc])); % Get local linearization
%     Bloc = double(subs(Ju,[x;u],[x_loc;u_loc]));
%     cloc = double(subs(f_ud_sym,[x;u],[x_loc;u_loc])) - Aloc*x_loc - Bloc*u_loc;
%     [U_loc,~,optval] = solveMPCprob(Aloc,Bloc,Cy,cloc,Q,R,Q,Np,-1, 1,[nan;nan;ymin;nan],[nan;nan;ymax;nan],x_loc,yr); % Get control input
%     u_loc = U_loc(1:m,1);
%     if(optval == Inf) % Detect infeasibility
%         ind_inf = [ind_inf i];
%         wasinfeas = 1;
%     end
%     x_loc = f_ud(0,x_loc,u_loc,deltaT,size(x_loc)); % Update true state
    
    % Store values
    XX_koop = [XX_koop x_koop];
    UU_koop = [UU_koop u_koop];
%     XX_loc = [XX_loc x_loc];
%     UU_loc = [UU_loc u_loc];
end

if(isempty(ind_inf))
    ind_inf = Nsim;
end

%% Plot (feedback control)
% Control signal
figure
p3 = plot([0:Nsim]*deltaT,ones(Nsim+1,1),'-k','linewidth',lw_koop-1); hold on
p4 = plot([0:Nsim]*deltaT,-ones(Nsim+1,1),'-k','linewidth',lw_koop-1); hold on
% p1 = plot([0:ind_inf(1)-1]*deltaT,UU_loc(1:ind_inf(1)),'--g','linewidth',lw_koop); hold on
p2 = plot([0:Nsim-1]*deltaT,UU_koop,'-b','linewidth',lw_koop); hold on
% axis([0,Tmax,min(UU_loc) - abs(min(UU_loc))*0.1, max(UU_loc)*1.1] )
LEG  = legend([p2,p3],'K-MPC','Constraint');
set(LEG,'Interpreter','latex','location','southeast')
set(gca,'FontSize',31);
 
figure
% Output (y = x_2)
p1=plot([0:Nsim]*deltaT,ymax*ones(Nsim+1,1),'-k','linewidth',lw_koop-1); hold on
p2=plot([0:Nsim]*deltaT,ymin*ones(Nsim+1,1),'-k','linewidth',lw_koop-1);
p3=plot([0:Nsim]*deltaT,XX_koop(1,:),'-b','linewidth',lw_koop); hold on
% p4=plot([0:ind_inf(1)-1]*deltaT,XX_loc(2,1:ind_inf(1)),'--g','linewidth',lw_koop);
p5=plot([0:Nsim-1]*deltaT, yrliftstack(1,:),'--r','linewidth',lw_koop);
LEG  = legend([p3,p5,p2],'K-MPC','Reference','Constraints');
set(LEG,'Interpreter','latex','location','southeast')
set(LEG,'Fontsize',18)
set(gca,'FontSize',20);

%%
figure
idx = 1;

plot([0:Nsim-1]*deltaT, XX_koop(idx,1:Nsim),'-b','linewidth',lw_koop); hold on
plot([0:Nsim-1]*deltaT, X(idx,1:Nsim),'-r','linewidth',lw_koop-1);
%%
function dYdt = pend(t,y,u)
    g = 9.81; 
    L =  0.365; 

    m = 0.0905; %mass of bob (kg)
    M = 1.12;  % mass of cart (kg)
    
    mux = 6.65;
    G = 7.5;

    x_ddot = G*u - mux*y(2,:) - m*L*y(4,:)*y(4,:)*sin(y(3,:));
    x_ddot = x_ddot / ( M+m*sin(y(3))* sin(y(3)) );

    theta_ddot =  G*u - mux*y(2,:) - m*L*y(4)*y(4,:)*sin(y(3,:))*cos(y(3,:)) + (M+m)*g*sin(y(3,:));
    theta_ddot = theta_ddot/ (L - m*L*cos(y(3,:))*cos(y(3,:)));
  
    dYdt = [ y(2,:); x_ddot; y(4,:); theta_ddot ];
end