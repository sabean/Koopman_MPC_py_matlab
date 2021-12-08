clear all
clc
close all

addpath('C:/Users/sabin/Desktop/KoopmanMPC-master/KoopmanMPC-master/Resources')
addpath('C:/Users/sabin/Desktop/KoopmanMPC-master/KoopmanMPC-master/Resources/qpOASES-3.1.0/interfaces/matlab') 


%% ****************************** Dynamics ********************************
%Runge-Kutta 4


n = 4; % Number of states
m = 1; % Number of control inputs
f_ud_pend_1 = @f_ud_pend_new;
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


%% Collect data
rng(115123)
disp('Starting data collection')
Ntraj = 1; Nsim = 1000;

Cy = [1 0 0 0; 0 0 1 0]; % Output matrix: y = Cy*x
ny = size(Cy,1); % Number of outputs

[X,Y,U] = CollectDataPendDelay(deltaT,uconstraints,stateconstrainsmin, stateconstrainsmax,Ntraj,Nsim,n);
fprintf('Data collection DONE \n');


%% Basis functions
basisFunction = 'rbf';
Nrbf = 100;
nD = 1;
n_zeta = (nD+1)*ny + nD*m; % dimension of the delay-embedded "state"
totalcols = size(X,2);
idx = randperm(totalcols);
cent = X(:,idx(1:Nrbf));

rbf_type = 'gauss';
liftFun = @(xx)( [xx;rbf(xx,cent,rbf_type)] );
%liftFun = @(x)( [ x ; ones(1,size(x,2)); -ones(1,size(x,2)) ; x(1:end-1,:).*(x(2:end,:)) ; x(1,:).*x(end,:) ; x.*x ]  ) ;

%Nlift = numel(liftFun(rand(n,1)));
Nlift =Nrbf + n_zeta;
%% Lift
disp('Starting LIFTING')

Xlift = liftFun(X);
Ylift = liftFun(Y);

% Regression

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
%fprintf( 'Regression residual %f \n', norm(Ylift - Alift*Xlift - Blift*U,'fro') / norm(Ylift,'fro') );


%% ************************* Predictor comparison *************************
Tmax = 5;
Nsim = Tmax/deltaT;

uprbs = ((umax-umin))*myprbs(Nsim,0.5)+umin;
u_dt = @(i)(  uprbs(i+1) );
%u_dt = @(i)((-1).^(round(i/30)));
f_cont_d = @(t,xx)( f_ud_pend_1(t,xx,u_dt(t), deltaT));


x0 = [0;0;pi;0];
x = x0;

% Delayed initial condition (assume random control input in the past)
xstart = [Cy*x ; NaN(nD*(ny+m),1)];
for i = 1:nD
    urand = 2*rand(m,1) - 1;
    xp = f_ud_pend_1(0,x,urand, deltaT);
    xstart = [Cy*xp ; urand; xstart(1:end-ny-m)];
    x = xp;
end


% Inital conditions
x_true = xp;
xlift = liftFun(xstart);

% Simulation
for i = 0:Nsim-1
    
    % True dynamics
    x_true = [x_true, f_ud_pend_1(0,x_true(:,end),u_dt(i),deltaT) ];
    
    % Koopman predictor
    xlift = [xlift Alift*xlift(:,end) + Blift*u_dt(i)];

end

%%
% figure
% stairs([0:Nsim-1]*deltaT,u_dt(0:Nsim-1),'linewidth',2); hold on
% title('Control input'); xlabel('time [s]')


lw_koop = 3;
s11= Cy*x_true;
sc11= Clift(1,:)*xlift;
sc12 = Clift(2,:)*xlift;


figure

plot([0:Nsim]*deltaT,s11(1,:),'--b','linewidth', lw_koop); hold on
plot([0:Nsim]*deltaT,sc11, '--r','linewidth',lw_koop)

LEG = legend('True','Koopman');
set(LEG,'Interpreter','latex','location','northeast','fontsize',30)
set(gca,'FontSize',25);



figure

plot([0:Nsim]*deltaT,s11(2,:),'--b','linewidth', lw_koop); hold on
plot([0:Nsim]*deltaT,sc12, '--r','linewidth',lw_koop)


LEG = legend('True','Koopman');
set(LEG,'Interpreter','latex','location','northeast','fontsize',30)
set(gca,'FontSize',25);


%% ********************** Feedback control ********************************
Tmax = 10 ;% Simlation legth
Nsim = Tmax/deltaT;

% Define Koopman controller
C = zeros(2,Nlift); C(2,3) = 1;% C(1,1) = 1;

% Weight matrices
Q = 120;
R = 0.001;
% Prediction horizon
Tpred = 0.1;
Np = round(Tpred / deltaT);
% Constraints
xlift_min = nan(Nlift,1);
xlift_max = nan(Nlift,1);
% Build Koopman MPC controller
koopmanMPC  = getMPC(Alift,Blift,C,0,Q,R,Q,Np,umin, umax, xlift_min, xlift_max,'qpoases');

REF = 'ref'; % 'step' or 'cos'
switch REF
    case 'step'
        ymin = -0.6;
        ymax = 0.6;
        yrr11 = 0.3*( -1 + 2*([1:Nsim] > Nsim/2)  ); % reference
        yrr12 = 0.3*( -1 + 2*([1:Nsim] > Nsim/2)  ); % reference
        
    case 'cos'
        ymin = -4;ymax = 4;
        ymin1 = -4;ymax1 =4 ;
        yrr11 = cos(2*pi*[1:Nsim] / Nsim); % reference
        yrr12 = cos(2*pi*[1:Nsim] / Nsim); % reference *(ymax-ymin)+ymin
    case 'ref'
        yrr11 = X(1,1:Nsim);
        yrr12 = X(2,1:Nsim); 
end
yrrc = [yrr11;yrr11;yrr12;yrr12];
x0=yrrc(:,1);
% Initial condition for the delay-embedded state (assuming zero control in the past)
x = x0;
zeta0 = [Cy*x ; NaN(nD*(ny+m),1)];
for i = 1:nD
    upast = zeros(m,1);
    xp = f_ud_pend_1(0,x,upast, deltaT);
    zeta0 = [Cy*xp ; upast ; zeta0(1:end-ny-m)]; 
    x = xp;
end
x0 = x;

x_koop = x0; 
zeta = zeta0; % Delay-embedded "state"

XX_koop = x0; UU_koop = [];

wasinfeas= 0;
ind_inf = [];

% Closed-loop simultion start
UU_koop1 = [];
XX_kp = [];
yrliftstack = [];
disp('Press any key for feedback control')
pause


tic
for i = 0:Nsim-1
    if(mod(i,10) == 0)
        fprintf('Closed-loop simulation: iterate %i out of %i \n', i, Nsim)
    end
    % Current value of the reference signal
    
    
    yr = yrrc(:,i+1);
    currhoz = (i+1)*deltaT;
   
    
    % Koopman MPC
    xlift = liftFun(zeta); % Lift
    yrrlift = [yr(1); yr(3)];
    yrliftstack = [yrliftstack yrrlift];
    u_koop = koopmanMPC(xlift,[yr(1);yr(3)]); % Get control input
    x_koop = f_ud_pend_1(0,x_koop,u_koop,deltaT); % Update true state
    
%     if(mod(currhoz, 1) == 0)
%         x_koop = yrrc(:,i+1); % Update true state 
%     else
%         x_koop = f_ud_pend(0,x_koop,u_koop,deltaT); % Update true state
%     end
    zeta = [ Cy*x_koop ; u_koop; zeta(1:end-ny-m)]; % Update delay-embedded state
    
    % Store values
    XX_koop = [XX_koop x_koop];
    UU_koop = [UU_koop u_koop];


end
toc

if(isempty(ind_inf))
    ind_inf = Nsim;
end

%% Plot (feedback control)
% Control signal
figure
p3 = plot([0:Nsim]*deltaT,ones(Nsim+1,1)*umin,'-k','linewidth',lw_koop-1); hold on
p4 = plot([0:Nsim]*deltaT,-ones(Nsim+1,1)*umax,'-k','linewidth',lw_koop-1); hold on
p2 = plot([0:Nsim-1]*deltaT,UU_koop,'-b','linewidth',lw_koop); hold on
LEG  = legend([p2,p3],'K-MPC','Constraint');
set(LEG,'Interpreter','latex','location','southeast')
set(gca,'FontSize',31);
%%
figure

p5=plot([0:Nsim-1]*deltaT, yrliftstack(1,:),'--r','linewidth',lw_koop);hold on
p3=plot([0:Nsim]*deltaT,XX_koop(1,:),'-b','linewidth',lw_koop);
LEG  = legend([p5,p3],'Reference','K-MPC');
set(LEG,'Interpreter','latex','location','southeast')
set(LEG,'Fontsize',18)
set(gca,'FontSize',20);
%%
figure

p5=plot([0:Nsim-1]*deltaT, yrliftstack(2,:),'--r','linewidth',lw_koop);hold on
p3=plot([0:Nsim]*deltaT,XX_koop(3,:),'-b','linewidth',lw_koop);
LEG  = legend([p5,p3],'Reference','K-MPC');
set(LEG,'Interpreter','latex','location','southeast')
set(LEG,'Fontsize',18)
set(gca,'FontSize',20);

