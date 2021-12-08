clear all
clc
close all

addpath('C:/Users/sabin/Desktop/KoopmanMPC-master/KoopmanMPC-master/Resources')
addpath('C:/Users/sabin/Desktop/KoopmanMPC-master/KoopmanMPC-master/Resources/qpOASES-3.1.0/interfaces/matlab') 


%% ****************************** Dynamics ********************************

n = 4; % Number of states
m = 1; % Number of control inputs

f_ud_pend = @f_ud_pend_new;

umin = -1;
umax = 1;


x1min = 0;
x1max = 0.5;

xdotmin = -1;
xdotmax = 1;

thetamin = -pi;
thetamax = pi;

thetadotmin = -10;
thetadotmax = 10;

uconstraints =[umin;umax];

stateconstrainsmin = [x1min;xdotmin;thetamin;thetadotmin];
stateconstrainsmax = [x1max;xdotmax;thetamax;thetadotmax];

% Discretize
deltaT = 0.01;


%% Collect data
rng(115123)
disp('Starting data collection')
Ntraj = 10; Nsim = 2000;

Cy = [1 0 0 0; 0 0 1 0]; % Output matrix: y = Cy*x
ny = size(Cy,1); % Number of outputs

[X,Y,U] = CollectDataPendTrig(deltaT,uconstraints,stateconstrainsmin, stateconstrainsmax,Ntraj,Nsim,n);
fprintf('Data collection DONE \n');

%% Basis functions
rng(115123)
Nrbf = 100;
RM1=[];RM2=[];RME1=[];RME2=[];
totalcols = size(X,2);
idx = randperm(totalcols);

pis = 2*pi*rand(100,1)-pi;
%rbfs = ["polyharmonic", "gauss", "invquad" ,"invmultquad", "thinplate"];
rbfs = [5, 10, 25, 50, 75, 100];
rbf_type = "invquad";
for kp = 1:length(rbfs)
    %rbf_type = rbfs(kp)
    
    RMSE1 =0;RMSE2=0;RMSEE1=0;RMSEE2=0;
    Nrbf = rbfs(kp);
    cent = X(:,idx(1:Nrbf));
    
   
    liftFun = @(xx)( [xx;rbf(xx,cent,rbf_type)] );

    Nlift = numel(liftFun(rand(n,1)));
    % Lift
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
    fprintf( 'Regression residual %f \n', norm(Ylift - Alift*Xlift - Blift*U,'fro') / norm(Ylift,'fro') );

    Tmax = 3;
    Nsim = Tmax/deltaT;

    uprbs = ((umax-umin))*myprbs(Nsim,0.5)+umin;
    u_dt = @(i)(  uprbs(i+1) );
   for j = 1:100

        
%        ************************* Predictor comparison *************************
        
        u_dt = @(i)((-1).^(round(i/30)));
        f_cont_d = @(t,xx)( f_ud_pend(t,xx,u_dt(t), deltaT));

        x0 = [0;0; pis(j);0];

        % Inital conditions
        %x0 = [0;0;pi;0];
        x_true = x0;
        xlift = liftFun(x0);


        % Simulation
        for i = 0:Nsim-1

            % True dynamics
            x_true = [x_true, f_ud_pend(0,x_true(:,end),u_dt(i),deltaT) ];

            % Koopman predictor
            xlift = [xlift Alift*xlift(:,end) + Blift*u_dt(i)];

        end
        RMSE1 = RMSE1 + sqrt(mean((x_true(1,:) - xlift(1,:)).^2));
        RMSE2 = RMSE2 + sqrt(mean((x_true(3,:) - xlift(3,:)).^2));
        RMSEE1 = RMSEE1 + 100* sqrt(mean((x_true(1,:) - xlift(1,:)).^2))/sqrt(mean((x_true(1,:).^2)));
        RMSEE2 = RMSEE2 + 100* sqrt(mean((x_true(3,:) - xlift(3,:)).^2))/sqrt(mean((x_true(3,:).^2)));
    end
    RM1 = [RM1 RMSE1/100]
    RM2 = [RM2 RMSE2/100]
    RME1 = [RME1 RMSEE1/100]
    RME2 = [RME2 RMSEE2/100]
end

%%

figure
stairs([0:Nsim-1]*deltaT,u_dt(0:Nsim-1),'linewidth',2); hold on
title('Control input'); xlabel('time [s]')


lw_koop = 3;
s11= Cy*x_true;
sc11= Clift(1,:)*xlift;
sc12 = Clift(3,:)*xlift;


figure

plot([0:Nsim]*deltaT,s11(1,:),'--b','linewidth', lw_koop); hold on
plot([0:Nsim]*deltaT,sc11, '--r','linewidth',lw_koop)

LEG = legend('True','Koopman');
set(LEG,'Interpreter','latex','location','northeast','fontsize',30)
set(gca,'FontSize',25);
%axis([0 1 -1.3 0.5])


figure

plot([0:Nsim]*deltaT,s11(2,:),'--b','linewidth', lw_koop); hold on
plot([0:Nsim]*deltaT,sc12, '--r','linewidth',lw_koop)


LEG = legend('True','Koopman');
set(LEG,'Interpreter','latex','location','northeast','fontsize',30)
set(gca,'FontSize',25);
%axis([0 1 -1.3 0.5])


%% ********************** Feedback control ********************************
Tmax = 10; % Simlation legth
Nsim = Tmax/deltaT;

x0 = [ 0; 0.; -pi; 0. ];

% Define Koopman controller
C = zeros(2,Nlift);  C(1,1) = 1;%C(2,3) = 1;

% Weight matrices
Q = 200;
R = 0.001;

% Prediction horizon
Tpred = 1;
Np = round(Tpred / deltaT);
% Constraints
xlift_min = [-0.8;nan(Nlift-1,1)];
xlift_max = [0.8;nan(Nlift-1,1)];

% xlift_min = nan(Nlift,1);
% xlift_max = nan(Nlift,1);

% Build Koopman MPC controller
koopmanMPC  = getMPC(Alift,Blift,C,0,Q,R,Q,Np,umin, umax, xlift_min, xlift_max,'qpoases');

yrliftstack = [];
REF = 'ref'; % 'step' or 'cos'
switch REF
    case 'step'
        ymin = -0.6;
        ymax = 0.6;
        yrr11 = 1*( -1 + 2*([1:Nsim] > Nsim/2)  ); % reference
        yrr12 = 1*( -1 + 2*([1:Nsim] > Nsim/2)  ); % reference
        
    case 'cos'
       ymin = -1;
        ymax = 1;
        yrr11 =0.8* cos(2*pi*[1:Nsim] / Nsim); % reference
        yrr12 =2*cos(2*pi*[1:Nsim] / Nsim); % reference
    case 'ref'
        yrr11 = X(1,1:Nsim);
        yrr12 = X(3,1:Nsim); 
end
yrrc = [yrr11;yrr11;yrr12;yrr12];

% Initial condition for the delay-embedded state (assuming zero control in the past)
x0=yrrc(:,1);
x = x0;

x_koop = x0; 
zeta = x0; % Delay-embedded "state"
XX_koop = x0; UU_koop = [];


% Closed-loop simultion start

disp('Press any key for feedback control')
pause
VIDEO =1;

for i = 0:Nsim-1
    if(mod(i,10) == 0)
        fprintf('Closed-loop simulation: iterate %i out of %i \n', i, Nsim)
    end
    yr = yrrc(:,i+1);
    currhoz = (i+1)*deltaT;
    

    
    % Current value of the reference signal
    %yr = yrr(i+1);
    
    % Koopman MPC
    xlift = liftFun(zeta); % Lift
    
    yrrlift = [yr(1); yr(3)];
    yrliftstack = [yrliftstack yrrlift];
    u_koop = koopmanMPC(xlift,[yr(1);yr(3)]); % Get control input
%     if(mod(currhoz, 1) == 0)
%         x_koop = yrrc(:,i+1); % Update true state 
%     else
%         x_koop = f_ud_pend(0,x_koop,u_koop,deltaT); % Update true state
%     end

    x_koop = f_ud_pend(0,x_koop,u_koop,deltaT); % Update true state
    zeta = x_koop; % Update delay-embedded state
    
    
    % Store values
    XX_koop = [XX_koop x_koop];
    UU_koop = [UU_koop u_koop];
end
%%
liftcxmin = 0.5;
liftcxmax = -0.5;
liftcthmin = 1.2;
liftcthmax = 4;

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
% Output (y = x_2)
%p1=plot([0:Nsim]*deltaT,liftcthmin*ones(Nsim+1,1),'-k','linewidth',lw_koop-1); hold on
%p2=plot([0:Nsim]*deltaT,liftcthmax*ones(Nsim+1,1),'-k','linewidth',lw_koop-1);
p3=plot([0:Nsim]*deltaT,XX_koop(3,:),'-b','linewidth',lw_koop); hold on
p5=plot([0:Nsim-1]*deltaT, yrliftstack(2,:),'--r','linewidth',lw_koop);
%LEG  = legend([p3,p5,p2],'K-MPC','Reference','Constraints');
LEG  = legend([p3,p5],'K-MPC','Reference');
set(LEG,'Interpreter','latex','location','southeast')
set(LEG,'Fontsize',18)
%axis([0,Tmax,liftcthmin,liftcthmax])
set(gca,'FontSize',20);
%%
figure
% Output (y = x_2)
%p1=plot([0:Nsim]*deltaT,liftcxmin*ones(Nsim+1,1),'-k','linewidth',lw_koop-1); hold on
%p2=plot([0:Nsim]*deltaT,liftcxmax*ones(Nsim+1,1),'-k','linewidth',lw_koop-1);
p3=plot([0:Nsim]*deltaT,XX_koop(1,:),'-b','linewidth',lw_koop); hold on
p5=plot([0:Nsim-1]*deltaT, yrliftstack(1,:),'--r','linewidth',lw_koop);
%LEG  = legend([p3,p5,p2],'K-MPC','Reference','Constraints');
LEG  = legend([p3,p5],'K-MPC','Reference');
set(LEG,'Interpreter','latex','location','southeast')
set(LEG,'Fontsize',18)
%axis([0,Tmax,liftcxmin,liftcxmax])
set(gca,'FontSize',20);

%%
%%save('PendulumTrajectoryData1','X','Y','U','Nsim','Ntraj','Xlift','Ylift', 'Clift')
v11=sqrt(mean((yrliftstack(1,1:999) - XX_koop(1,1:999)).^2))

v21=100*norm(yrliftstack(1,1:999) - XX_koop(1,1:999))/norm(yrliftstack(1,1:999))
