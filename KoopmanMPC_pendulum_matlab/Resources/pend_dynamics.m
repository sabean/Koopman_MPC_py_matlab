function dYdt = pend_dynamics(t,y,u)
    g = 9.81; 
    L =  0.365; 

    m = 0.0905; %mass of bob (kg)
    M = 1.12;  % mass of cart (kg)
    
    mux = 6.65;
    G = 7.5;

    x_ddot = G*u - mux*y(2,:) - m*L*y(4,:)*y(4,:)*sin(y(3,:));
    x_ddot = x_ddot / ( M+m*sin(y(3,:))* sin(y(3,:)) );

    theta_ddot =  G*u - mux*y(2,:) - m*L*y(4,:)*y(4,:)*sin(y(3,:))*cos(y(3,:)) + (M+m)*g*sin(y(3,:));
    theta_ddot = theta_ddot/ (L - m*L*cos(y(3,:))*cos(y(3,:)));
  
    dYdt = [ y(2,:); x_ddot; y(4,:); theta_ddot ];
end