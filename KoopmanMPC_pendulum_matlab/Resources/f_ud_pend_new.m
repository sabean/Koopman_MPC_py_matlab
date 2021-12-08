function yout = f_ud_pend_new(t, y0, u, dt)
    tspan = [0, dt];
    [t,y] = ode45(@pend_dynamics, tspan, y0, [], u);
    yout = y(end,:).';
end

function dYdt = pend_dynamics(t,y,u)
    g = 9.81; 
    L =  0.365; 

    m = 0.0905; %mass of bob (kg)
    M = 1.12;  % mass of cart (kg)
    
    mux = 6.65;
    G = 7.5;

    %x_ddot = G*u - mux*y(2) - m*L*y(4)*y(4)*sin(y(3)) - m*g*cos(y(3))*sin(y(3));
    x_ddot = G*u - mux*y(2) - m*L*y(4)*y(4)*sin(y(3));
    x_ddot = x_ddot / ( M+m*sin(y(3))* sin(y(3)) );

    theta_ddot =  G*u - mux*y(2) - m*L*y(4)*y(4)*sin(y(3))*cos(y(3)) + (M+m)*g*sin(y(3));
    theta_ddot = theta_ddot/ (L - m*L*cos(y(3))*cos(y(3)));
    %theta_ddot =  -cos(y(3))*G*u + mux*y(2)*cos(y(3)) - m*L*y(4)*y(4)*sin(y(3))*cos(y(3)) + (M+m)*g*sin(y(3));
    %theta_ddot = theta_ddot/ (L*(M + m*sin(y(3))*sin(y(3))));
    dYdt = [ y(2); x_ddot; y(4); theta_ddot ];
end