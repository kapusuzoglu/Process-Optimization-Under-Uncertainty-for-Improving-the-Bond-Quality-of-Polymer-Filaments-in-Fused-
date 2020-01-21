function[bond,t_bond,theta]=Neck_Growth_Model(T,dt,tfinal,w,h)

% ======================================================================= %
% automatic step size Runge-Kutta-Fehlberg integration method (Burden and Faires, 1985)
% To overcome numerical instabilities when theta = 0, the initial boundary condition
% is fixed at a time value slightly different than zero and the corresponding 
% value of theta is determined from Eq. 15

% They found that the majority of neck growth and sintering in the ABS fibers
% occurred when the interphase temperature was above 200°C, which meant (based on heat
% transfer analysis and confirmed by experiments) that the nozzle temperature had a large
% effect on sintering while environment temperature had little effect.

% In which delta is time step. For this case, delta is set equal to 2*dt. 
% dt is the time step that is used for the interval loop
% ======================================================================= %


% Get Temperature and time for the cooling
% [T t]=cooling_model_ABS(To,Tinf,dt,tfinal);

% ao = w/2; % initial radius: half of layer width in meters
ao = sqrt(w/2*h/2); % initial radius: half of layer width in meters

%Surface tension
% gamma=.0342; % polycarbonate
gamma=.029; % N/m for ABS P400 at 240 celcius
% gamma=.047; % N/m for ABS P400 at 220 celcius
% gamma=.054; % N/m for ABS P400 at 200 celcius
% with a temp. dependent of Delta Gamma/ Delta T = - 0.000345 N/m/K
Delta_Gamma = -0.00345*ones(1,ceil(tfinal/dt)+1);
KK = 273.15;
[~,idx]=min(abs(T-(240+KK)));
if T(1)>240+KK
    Delta_Gamma(1:idx) = -0.0005;
% Delta_Gamma(1:idx) = -0.000000001;
end

% Temperature dependent viscosity
eta_r=5100; % %Viscosity at temp 240 celc
T_r_C = 240; % reference temperature in C
T_r = T_r_C+KK; % in K
b_m = 0.056; % model parameter for temp dependent viscosity

Eta(1) = eta_r*exp(-b_m*(T(1)-T_r));
Gamma(1) = gamma + Delta_Gamma(1) * (T(1)-T_r);

Eta(2) = eta_r*exp(-b_m*(T(2)-T_r));
Gamma(2) = gamma + Delta_Gamma(2) * (T(2)-T_r);

theta(1)=0;
t_bond(1)=0;

theta(2)=sqrt(2*(dt)*Gamma(2)/(Eta(2)*ao)); % Eq. 15
t_bond(2)=2*dt;

for j=4:2:(tfinal/dt-1)
    delta_t=2*dt;
    t_bond(j/2+1)=t_bond(j/2)+delta_t;
    %k1 calculation at t_bond(i/2)
    eta_1=eta_r*exp(-b_m*(T(j-1)-T_r));
    gamma_1 = gamma+ Delta_Gamma(j-1) * (T(j-1)-T_r);
    theta_1=theta(j/2);
    k1=(gamma_1/(3*ao*eta_1*sqrt(pi)))*((pi-theta_1)*cos(theta_1)+...
        sin(theta_1))*((pi-theta_1+sin(theta_1)*(cos(theta_1)))^(1/2))/...
        (((pi-theta_1)^2)*((sin(theta_1))^2));
    
    %k2 calculation
    eta_2=eta_r*exp(-b_m*(T(j)-T_r));
    gamma_2 = gamma+ Delta_Gamma(j) * (T(j)-T_r);
    theta_2=theta(j/2)+dt*k1;
    k2=(gamma_2/(3*ao*eta_2*sqrt(pi)))*((pi-theta_2)*cos(theta_2)+...
        sin(theta_2))*((pi-theta_2+sin(theta_2)*(cos(theta_2)))^(1/2))/...
        (((pi-theta_2)^2)*((sin(theta_2))^2));
    %k3 calculation
    eta_3=eta_2;
    gamma_3=gamma_2;
    theta_3=theta(j/2)+dt*k2;
    k3=(gamma_3/(3*ao*eta_3*sqrt(pi)))*((pi-theta_3)*cos(theta_3)+...
        sin(theta_3))*((pi-theta_3+sin(theta_3)*(cos(theta_3)))^(1/2))/...
        (((pi-theta_3)^2)*((sin(theta_3))^2));
    %k4 calculation
    eta_4=eta_r*exp(-b_m*(T(j+1)-T_r));
    gamma_4 = gamma+ Delta_Gamma(j+1) * (T(j+1)-T_r);
    theta_4=theta(j/2)+2*dt*k3;
    k4=(gamma_4/(3*ao*eta_4*sqrt(pi)))*((pi-theta_4)*cos(theta_4)+...
        sin(theta_4))*((pi-theta_4+sin(theta_4)*(cos(theta_4)))^(1/2))/...
        (((pi-theta_4)^2)*((sin(theta_4))^2));
    
    % Temperature dependent viscosity and surface tension
    Eta(j/2+1) = eta_2;
    Gamma(j/2+1) = gamma_2;
    
    % theta
    theta(j/2+1)=theta(j/2)+(1/6)*delta_t*(k1+2*k2+2*k3+k4);
end
% y = a*sin(theta)
bond=ao*sin(theta);

% % At time t with radius r
% r = ao*sqrt(pi)./(sqrt(pi-theta+sin(theta).*cos(theta)));
% 
% bond=r.*sin(theta);


% Gamma
% Delta_Gamma
% dimensionless time Dim_t
% Dim_t = t_bond .* Gamma./(Eta*ao);

% Todd = T(1:2:end-1);          % Odd-Indexed Elements




% figure
% plot(Dim_t,sin(theta),'-k','linewidth',2);
% xlabel('Dimensionless time $(t\Gamma/\eta a_0)$', 'FontName', 'Times New Roman', ...
%                     'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
% ylabel('Dimensionless neck radius $(y/a_0)$', 'FontName', 'Times New Roman', ...
%                     'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
% set(gcf, 'PaperPosition', [0 0 4 6]); %Position plot at left hand corner with width 5 and height 5.
% set(gcf, 'PaperSize', [4 6]); %Set the paper to have width 5 and height 5.










% figure
% yyaxis left
% plot(t_bond,bond*1e3,'-b','linewidth',2);hold on
% xlabel('Time (s)', 'FontName', 'Times New Roman', ...
%                     'FontSize',16,'Color','k', 'Interpreter', 'tex')
% ylabel('Neck radius (mm)', 'FontName', 'Times New Roman', ...
%                     'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
% yyaxis right
% plot(t_bond,Todd-KK,'linewidth',2);
% ylabel('Temperature ({\circ}C)', 'FontName', 'Times New Roman', ...
%                     'FontSize',16,'Color','k', 'Interpreter', 'tex')
% set(gcf, 'PaperPosition', [0 0 4 6]); %Position plot at left hand corner with width 5 and height 5.
% set(gcf, 'PaperSize', [4 6]); %Set the paper to have width 5 and height 5.
% legend('Neck radius predictions','Interface temperature')




% figure
% plot(Todd-KK,bond*1e3,'-b','linewidth',2);hold on
% xlabel('Temperature ({\circ}C)', 'FontName', 'Times New Roman', ...
%                     'FontSize',16,'Color','k', 'Interpreter', 'tex')
% ylabel('Neck radius (mm)', 'FontName', 'Times New Roman', ...
%                     'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
% set(gcf, 'PaperPosition', [0 0 4 6]); %Position plot at left hand corner with width 5 and height 5.
% set(gcf, 'PaperSize', [4 6]); %Set the paper to have width 5 and height 5.
% h = gca;  % Handle to currently active axes
% set(h, 'XDir', 'reverse');


% figure
% % eta_r=5100; % temp 240 celc
% % T_r=513;
% % T = linspace(540,370,10);
% % b_m =0.056; % model parameter
% % Eta=eta_r.*exp(-b_m*(T-T_r));
% plot(Todd,Eta)

% fname = 'C:\Users\berkc\Dropbox\Vandy\Research\AM\MatlabCode\Figs';
% set(gca(), 'LooseInset', get(gca(), 'TightInset'));
% saveas(gcf, fullfile(fname,'\DeterProg_RDO'), 'pdf') %Save figure 

end