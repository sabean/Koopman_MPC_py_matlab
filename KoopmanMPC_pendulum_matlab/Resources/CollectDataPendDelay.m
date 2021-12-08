function [X,Y,U] = CollectDataPendDelay(dt,uconstraints,stateconstrainsmin, stateconstrainsmax,Ntraj,SimLength, n)
% Collect data for Pendulum equation
rng(115123)
disp('Starting data collection ...')

nD = 1; % Number of delays
Cy = [1 0 0 0; 0 0 1 0]; % Output matrix: y = Cy*x
ny = size(Cy,1); % Number of outputs
m = 1;
% random input

% Transition mapping of the controlled dynamical system
f = @(x,u)(f_ud_pend_new(0, x,u,dt));


% 

% initialize 
X = []; Y = []; U=[];

stt = 0;
ste = 10;
% loop over trajectories
for i = 1:Ntraj
    % Intial state is constrainsted
%     x11 = rand(1,1)* (x1max-x1min)+x1min;
%     x12 = rand(1,1)*(xdotmax-xdotmin)+xdotmin;
%     theta11 = rand(1,1)* (thetamax-thetamin)+thetamin;
%     theta12 = rand(1,1)*(thetadotmax-thetadotmin)+thetadotmin;
%     xcurr = [x11;x12;theta11;theta12];
    xcurr =[0;0;pi;0];
    xx = [];
    yy = [];
    tic
    fprintf('Trajectory %d out of %d \n',i,Ntraj)
    % loop over each time step
    t = linspace(stt, ste, SimLength);
    omega = (3*pi- pi)*rand(1)+pi;
    phi = (2*pi)*rand(1);
    AMP = (1- 0.1)*rand(1)+0.1;
    uu = AMP*sin(omega*t + phi);
    stt = stt+ste;
    ste = ste+ste;
    zeta_current = [Cy*xcurr; NaN(nD*(ny+m),1)];
        
    for j = 1:SimLength 
        xnext = f(xcurr(:,end),uu(j));
        zeta_prev = zeta_current;
        zeta_current = [[Cy*xnext ; uu(j)] ;  zeta_current( 1:end-ny-m , : ) ];
        if(j > nD)
            xx = [xx zeta_prev];
            yy = [yy zeta_current];
            U  = [U,uu(j)];
        end
        xcurr = xnext;
        
        % if the solution diverges, go to the next trajectory
%         if ~isempty(find(isnan(xx(:,end)),1))
%             break
%         end
        
    end
    toc
    % Store
    X = [X xx(:,1:end)];
    Y = [Y yy(:,1:end)];

end


save('PendulumTrajectoryData','X','Y','U','SimLength','Ntraj')
end