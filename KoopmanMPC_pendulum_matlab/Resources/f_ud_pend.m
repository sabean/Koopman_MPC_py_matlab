function yout = f_ud_pend(t, y0, u, dt)
    tspan = [0, dt];
    [t,y] = ode45(@pend_dynamics, tspan, y0, [], u);
    yout = y(end,:).';
end

function dYdt = pend_dynamics(t,y,u)
    g = 9.8; 
    L = 0.22; 

    m = 0.41; %mass of bob (kg)
    M = 2.0;  % mass of cart (kg)


    x_ddot = u - m*L*y(4)*y(4) * cos( y(3) ) + m*g*cos(y(3)) *  sin(y(3));
    x_ddot = x_ddot / ( M+m-m* sin(y(3))* sin(y(3)) );

    theta_ddot = -g/L * cos( y(3) ) - 1./L * sin( y(3) ) * x_ddot;

    damping_theta =  0; %- 0.5*y[3]
    damping_x =  0; %- 1.0*y[1]
    dYdt = [ y(2); x_ddot + damping_x; y(4); theta_ddot + damping_theta ];
end