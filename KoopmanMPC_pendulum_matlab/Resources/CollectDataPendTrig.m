function [X,Y,U] = CollectDataPendTrig(dt,uconstraints,stateconstrainsmin, stateconstrainsmax,Ntraj,SimLength, n)
% Collect data for Pendulum equation
disp('Starting data collection ...')


% random input
umin = uconstraints(1); umax =  uconstraints(2);
% Ubig= rand(1,SimLength,Ntraj)* (umax-umin)+umin;
% 
% 
% period =1.5;

% Transition mapping of the controlled dynamical system
f = @(x,u)(f_ud_pend_new(0, x,u,dt));

% x1min = stateconstrainsmin(1);xdotmin = stateconstrainsmin(2);thetamin = stateconstrainsmin(3);thetadotmin = stateconstrainsmin(4);
% x1max = stateconstrainsmax(1);xdotmax = stateconstrainsmax(2);thetamax = stateconstrainsmax(3);thetadotmax = stateconstrainsmax(4);


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
    %xx = [x11;x12;theta11;theta12];
    xx =[0;0;pi;0];
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
    for j = 1:SimLength 
        xx = [xx f(xx(:,end),uu(j))];
        U  = [U,uu(j)];
        % if the solution diverges, go to the next trajectory
        if ~isempty(find(isnan(xx(:,end)),1))
            break
        end
        
    end
    toc
    % Store
    X = [X xx(:,1:end-1)];
    Y = [Y xx(:,2:end)];

end


save('PendulumTrajectoryData','X','Y','U','SimLength','Ntraj')
end