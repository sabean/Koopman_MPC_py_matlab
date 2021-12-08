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
f_u = @pend_dynamics;
k1 = @(t,x,u) (  f_u(t,x,u) );
k2 = @(t,x,u) ( f_u(t,x + k1(t,x,u)*deltaT/2,u) );
k3 = @(t,x,u) ( f_u(t,x + k2(t,x,u)*deltaT/2,u) );
k4 = @(t,x,u) ( f_u(t,x + k1(t,x,u)*deltaT,u) );
f_pp= @(t,x,u) ( x + (deltaT/6) * ( k1(t,x,u) + 2*k2(t,x,u) + 2*k3(t,x,u) + k4(t,x,u)  )   );


%% Collect data
rng(115123)
disp('Starting data collection')
Ntraj = 20; Nsim = 1000;

Cy = [1 0 0 0; 0 0 1 0]; % Output matrix: y = Cy*x
ny = size(Cy,1); % Number of outputs

[X,Y,U] = CollectDataPendDelay(deltaT,uconstraints,stateconstrainsmin, stateconstrainsmax,Ntraj,Nsim,n);
fprintf('Data collection DONE \n');

%% Basis functions
basisFunction = 'rbf';
RMSE = [];
selrbf = [5 10 25 50 75 100 150 200 500];
nD = 1;
n_zeta = (nD+1)*ny + nD*m; % dimension of the delay-embedded "state"
totalcols = size(X,2);
idx = randperm(totalcols);
rbf_type = 'polyharmonic';
pis = 2*pi*rand(100,1)-pi;
for kp = 1:length(selrbf)
    Nrbf = selrbf(kp);
    cent = X(:,idx(1:Nrbf));


    liftFun = @(xx)( [xx;rbf(xx,cent,rbf_type)] );
    Nlift =Nrbf + n_zeta;
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


    % ************************* Predictor comparison *************************
    Tmax = 5;
    Nsim = Tmax/deltaT;

    uprbs = ((umax-umin))*myprbs(Nsim,0.5)+umin;
    u_dt = @(i)(  uprbs(i+1) );
    %u_dt = @(i)((-1).^(round(i/30)));
    f_cont_d = @(t,xx)( f_ud_pend_1(t,xx,u_dt(t), deltaT));
    RMSE1 = 0;RMSE2 = 0;

    for j = 1:100

        %x0 = rand(4,1)-0.5;
        x0 = [0;0; pis(j);0];
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
        v1= Cy*x_true;
        v11 = v1(1,:);
        v12 = v1(2,:);

        v21= Clift(1,:)*xlift;
        v22 = Clift(2,:)*xlift;

        RMSE1 = RMSE1 + 100*norm(v21-v11)/norm(v11);
        RMSE2 = RMSE2 + 100*norm(v22-v12)/norm(v12);

    end
    RMSE1 = RMSE1/100;
    RMSE2 = RMSE2/100;
    RMSE =[RMSE (RMSE1+RMSE2)/2]
end 

%%

figure
stairs([0:Nsim-1]*deltaT,u_dt(0:Nsim-1),'linewidth',2); hold on
title('Control input'); xlabel('time [s]')


lw_koop = 3;
s11= Cy*x_true;
sc11= Clift(1,:)*xlift;
sc12 = Clift(2,:)*xlift;
%s13= Cy*xloc;

figure

plot([0:Nsim]*deltaT,s11(1,:),'--r','linewidth', lw_koop); hold on
plot([0:Nsim]*deltaT,sc11, '--b','linewidth',lw_koop)
%plot([0:Nsim]*deltaT,s13(1,:), '--c','linewidth',lw_koop-1)
LEG = legend('True','Koopman','Local at x0');
set(LEG,'Interpreter','latex','location','northeast','fontsize',30)
set(gca,'FontSize',25);


figure

plot([0:Nsim]*deltaT,s11(2,:),'-r','linewidth', lw_koop); hold on
plot([0:Nsim]*deltaT,sc12, '--b','linewidth',lw_koop)
%plot([0:Nsim]*deltaT,s13(2,:), '--c','linewidth',lw_koop-1)

LEG = legend('True','Koopman', 'Local at x0');
set(LEG,'Interpreter','latex','location','northeast','fontsize',30)
set(gca,'FontSize',25);




%% ********************** Feedback control ********************************
Tmax = 6; % Simlation legth
Nsim = Tmax/deltaT;


% Define Koopman controller
C = zeros(2,Nlift);  C(1,1) = 1;  C(2,3) = 1;

% Weight matrices
Q = 120;
R = 0.01;
% Prediction horizon
Tpred = 1;
Np = round(Tpred / deltaT);
% Constraints
xlift_min = [-5;nan(Nlift-1,1)];
xlift_max = [2;nan(Nlift-1,1)];

% xlift_min = liftFun(stateconstrainsmin);
% xlift_max = liftFun(stateconstrainsmax);
% Build Koopman MPC controller
koopmanMPC  = getMPC(Alift,Blift,C,0,Q,R,Q,Np,umin, umax, xlift_min, xlift_max,'qpoases');

REF = 'step'; % 'step' or 'cos'
switch REF
    case 'step'
        ymin = -0.6;
        ymax = 0.6;
        yrr11 = 0.3*( -1 + 2*([1:Nsim] > Nsim/2)  ); % reference
        yrr12 = 0.3*( -1 + 2*([1:Nsim] > Nsim/2)  ); % reference
        
    case 'cos'
        ymin = 2.6;ymax = 3.6;
        ymin1 = -0.8;ymax1 =0.8 ;
        yrr11 = cos(2*pi*[1:Nsim] / Nsim)*(ymax1-ymin1)+ymin1; % reference
        yrr12 = cos(2*pi*[1:Nsim] / Nsim)*(ymax-ymin)+ymin; % reference
    case 'ref'
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
    xp = f_ud_pend_1(0,x,upast, deltaT);
    zeta0 = [Cy*xp ; upast ; zeta0(1:end-ny-m)]; 
    x = xp;
end
x0 = x;

x_koop = x0; 
zeta = zeta0; % Delay-embedded "state"

XX_koop = x0; UU_koop = [];
% %%
% % Get Jacobian of the true dynamics (for local linearization MPC)
% g = 9.81; L =  0.365; m = 0.0905; M = 1.12; mux = 6.65; G = 7.5;
% 
% xddot = @(t,y,u)((G*u - mux*y(2) - m*L*y(4)*y(4)*sin(y(3)))/( M+m*sin(y(3))* sin(y(3)) ));
% theta_ddot = @(t,y,u)((G*u - mux*y(2) - m*L*y(4)*y(4)*sin(y(3))*cos(y(3)) + (M+m)*g*sin(y(3)))/(L - m*L*cos(y(3))*cos(y(3))));
% pend_dynamics = @(t,y,u)([ y(2); x_ddot(t,y,u); y(4); theta_ddot(t,y,u) ]);
% f_pend = @(t,y,u)(ode45(@pend_dynamics, [0 deltaT], y, [], u));
% f_res = @(y,u)(f_pend(0,y,u).y(:,end));

x_loc = x0;
XX_loc = x0; UU_loc = [];
u_loc = 0;
x = sym('x',[4 1]); syms u;
f_ud_sym = f_pp(0,x,u);
Jx = jacobian(f_ud_sym,x);
Ju = jacobian(f_ud_sym,u);


wasinfeas= 0;
ind_inf = [];

% Closed-loop simultion start
UU_koop1 = [];

yrliftstack = [];
disp('Press any key for feedback control')
pause
m=1;
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
    u_koop = koopmanMPC(xlift,[yr(1);yr(3)]); % Get control input
    x_koop = f_ud_pend_1(0,x_koop,u_koop,deltaT); % Update true state
    zeta = [ Cy*x_koop ; u_koop; zeta(1:end-ny-m)]; % Update delay-embedded state
    
    % Local linearization MPC
    Aloc = double(subs(Jx,[x;u],[x_loc;u_loc])); % Get local linearization
    Bloc = double(subs(Ju,[x;u],[x_loc;u_loc]));
    cloc = double(subs(f_ud_sym,[x;u],[x_loc;u_loc])) - Aloc*x_loc - Bloc*u_loc;
    [U_loc,~,optval] = solveMPCprob(Aloc,Bloc,Cy,cloc,Q,R,Q,Np,-1, 1,[-5;nan;nan;nan],[2;nan;nan;nan],x_loc,yrrlift); % Get control input
    u_loc = U_loc(1:m,1);
    if(optval == Inf) % Detect infeasibility
        ind_inf = [ind_inf i];
        wasinfeas = 1;
    end
  
    x_loc = f_pp(0,x_loc,u_loc); % Update true state
    
    
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
%p1=plot([0:Nsim]*deltaT,x1min*ones(Nsim+1,1),'-k','linewidth',lw_koop-1); hold on
%p2=plot([0:Nsim]*deltaT,x1max*ones(Nsim+1,1),'-k','linewidth',lw_koop-1);
p3=plot([0:Nsim]*deltaT,XX_koop(3,:),'-b','linewidth',lw_koop); hold on
p5=plot([0:Nsim-1]*deltaT, yrliftstack(2,:),'--r','linewidth',lw_koop);
p4 = plot([0:ind_inf(1)-1]*deltaT,XX_loc(3,1:ind_inf(1)),'--g','linewidth',lw_koop);
LEG  = legend([p3,p5,p4],'K-MPC','Reference','L-MPC');
set(LEG,'Interpreter','latex','location','southeast')
set(LEG,'Fontsize',18)
set(gca,'FontSize',20);