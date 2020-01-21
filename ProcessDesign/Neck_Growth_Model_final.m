
function[BondLengths] = Neck_Growth_Model(T,inputs,ProcessParam,Ref_temp,...
                        tfinal,dt,time,firstEl,b2_Samples)

    % Neck growth calculation using the experiment temperature
    % ==================================================================== %
    % automatic step size Runge-Kutta-Fehlberg integration method
    % (Burden and Faires, 1985)
    % To overcome numerical instabilities when theta = 0, the initial BC.
    % is fixed at a time value slightly different than zero and the
    % corresponding value of theta is determined from Eq. 15
    
    % They found that the majority of neck growth and sintering in the ABS
    % fibers occurred when the interphase temperature was above 200°C,
    % which meant (based on heat transfer analysis and confirmed by
    % experiments) that the nozzle temperature had a large effect on
    % sintering while environment temperature had little effect.
    
    % In which delta is time step. For this case, delta is set equal to 2*dt.
    % dt is the time step that is used for the interval loop
    % ==================================================================== %

    %Filament dimensions
    w = ProcessParam(1); % Layer Thickness (meters)
    h = ProcessParam(2); % Layer Height (meters)
%     layer = ProcessParam(3); % considered layer number
%     numFilaments = ProcessParam(4); % considered interface number on that layer

    % ao = w/2; % initial radius: half of layer width in meters
    %     ao = sqrt(w/2*h/2) % initial radius: half of layer width in meters
    aa=w/2; bb=h/2;
    % initial radius of an ellipse at 45 degrees
    ao=aa*bb/sqrt(aa^2*sin(0.25*pi)^2+bb^2*cos(0.25*pi)^2); % in [m]
    % Material properties - ABS
    %Surface tension
    gamma = inputs(1);
    % gamma=.029; % N/m for ABS P400 at 240 celcius
    % gamma=.047; % N/m for ABS P400 at 220 celcius
    % gamma=.054; % N/m for ABS P400 at 200 celcius
    % with a temp. dependent of Delta Gamma/ Delta T = - 0.000345 N/m/K
    b1=inputs(2);
    % b1 = 0.00345;
    Delta_Gamma = -b1*ones(1,ceil(tfinal/dt)+1);
    
    % Temperature dependent viscosity
    % eta_r=5100; % %Viscosity at temp 240 celc
    eta_r = inputs(3);
    % b2 = 0.056; % model parameter for temp dependent viscosity
    b2 = b2_Samples; % model parameter for temp dependent viscosity
    
    Kelv = 273.15; % Kelvin conversion
    [~,idx]=min(abs(T-(240+Kelv)));
    if T(1)>240+Kelv
        Delta_Gamma(1:idx) = -0.0005;
        % Delta_Gamma(1:idx) = -0.000000001;
    end
    
    T_r_C = 240; % reference temperature in C
    T_r = T_r_C+Kelv; % in K
    
    Eta(1) = eta_r*exp(-b2*(T(1)-T_r));
    Gamma(1) = gamma + Delta_Gamma(1) * (T(1)-T_r);
    
    Eta(2) = eta_r*exp(-b2*(T(2)-T_r));
    Gamma(2) = gamma + Delta_Gamma(2) * (T(2)-T_r);
    
    theta(1)=0;
    t_bond(1)=0;
    
    theta(2)=sqrt(2*(dt)*Gamma(2)/(Eta(2)*ao)); % Eq. 15
    t_bond(2)=2*dt;
    
    for jjj=4:2:(tfinal/dt-1)
        delta_t=2*dt;
        t_bond(jjj/2+1)=t_bond(jjj/2)+delta_t;
        %k1 calculation at t_bond(i/2)
        eta_1=eta_r*exp(-b2*(T(jjj-1)-T_r));
        gamma_1 = gamma+ Delta_Gamma(jjj-1) * (T(jjj-1)-T_r);
        theta_1=theta(jjj/2);
        k1=(gamma_1/(3*ao*eta_1*sqrt(pi)))*((pi-theta_1)*cos(theta_1)+...
            sin(theta_1))*((pi-theta_1+sin(theta_1)*(cos(theta_1)))^(1/2))/...
            (((pi-theta_1)^2)*((sin(theta_1))^2));
        
        %k2 calculation
        eta_2=eta_r*exp(-b2*(T(jjj)-T_r));
        gamma_2 = gamma+ Delta_Gamma(jjj) * (T(jjj)-T_r);
        theta_2=theta(jjj/2)+dt*k1;
        k2=(gamma_2/(3*ao*eta_2*sqrt(pi)))*((pi-theta_2)*cos(theta_2)+...
            sin(theta_2))*((pi-theta_2+sin(theta_2)*(cos(theta_2)))^(1/2))/...
            (((pi-theta_2)^2)*((sin(theta_2))^2));
        %k3 calculation
        eta_3=eta_2;
        gamma_3=gamma_2;
        theta_3=theta(jjj/2)+dt*k2;
        k3=(gamma_3/(3*ao*eta_3*sqrt(pi)))*((pi-theta_3)*cos(theta_3)+...
            sin(theta_3))*((pi-theta_3+sin(theta_3)*(cos(theta_3)))^(1/2))/...
            (((pi-theta_3)^2)*((sin(theta_3))^2));
        %k4 calculation
        eta_4=eta_r*exp(-b2*(T(jjj+1)-T_r));
        gamma_4 = gamma+ Delta_Gamma(jjj+1) * (T(jjj+1)-T_r);
        theta_4=theta(jjj/2)+2*dt*k3;
        k4=(gamma_4/(3*ao*eta_4*sqrt(pi)))*((pi-theta_4)*cos(theta_4)+...
            sin(theta_4))*((pi-theta_4+sin(theta_4)*(cos(theta_4)))^(1/2))/...
            (((pi-theta_4)^2)*((sin(theta_4))^2));
        
        % theta
        theta(jjj/2+1)=theta(jjj/2)+(1/6)*delta_t*(k1+2*k2+2*k3+k4);
    end
    % y = a*sin(theta)
    
    % bond length = initial radius*sin(angle)
    bondLength=ao*sin(theta);
    
    % find the time when temperature reaches 112 celcius
    index1=find(Ref_temp+Kelv>T,1,'first');
    
%     if rr==1
%         % find the time index that corresponds to the same time
%         [~,index2]=min(abs(time(index1)-t_bond));
%     else
%         % find the time index that corresponds to the same time
%         [~,index2]=min(abs(time(index1+firstEl-1)-t_bond));
%     end

    % find the time index that corresponds to the same time
    if firstEl~=0
        [~,index2]=min(abs(time(index1+firstEl-1)-t_bond));
    else
        [~,index2]=min(abs(time(index1)-t_bond));
    end
    
    % save neck radius for each interface on a given layer
    BondLengths = 2*bondLength(index2)*1e3; % neck radius in mm
    
end  