% function[AllBondLengths] = Experiment_FirstLayer(ProcessParam,inputs,NumTest,Ref_temp)
    
clc;clear;

format long

% Load Gaussian Process model for K1 and Work done
GPStruc_K1 = load('06 GP Model - K1 032 5-90\Fitting\gprMdl_K1');
GPStruc_W = load('06 GP Model - W 032 5-90\Fitting\gprMdl_W');
GP = {GPStruc_K1 GPStruc_W};


% ======================================================================= %
%__________________________       INPUTS     ______________________________
%Process Variables
%__________________________________________________________________________
%Filament dimensions
wth = 0.8e-3; % layer width [meters]
hth = 0.7e-3; % layer height [meters]
layer = 1; % which layer we are looking at for the bond lengths 
% numFilaments-1: (numFilaments-1)th interface on that given layer
numFilaments = 15;

% Process parameters
ProcessParam = [wth, hth, layer, numFilaments];

% 2 layers & 15 filaments
L = 0.035; % filament lengths in meters
Matrix = ones(2,15);
x = L/2; % location of the cut in meters



% Material parameters for the bond length model
Surf_tens = 0.29; % Surface tension
b1_ = 0.00345; % Model parameter for temp. dependent surface tension
Visco = 5100; % viscosity
% b2_ = 0.0193; % for temp above 130 celcius
b2_=0.0016;

Ref_Temp = 130; %celcius

inputs = [Surf_tens, b1_, Visco, b2_];


Kel = 273.15; % 0 Celcius in K
% convert Celcius to Kelvin

NumTest=1;

% ======================================================================= %
%% Experiment
%% This will plot the temperatures at line B (x = 175 pixels)
% Row pixels that are manually chosen which represents the interface 
% between lines and column pixels represent point B
row=[41,38,35,32,30,27,25,22,20,18,15,13,10,8];
% start_frame=19; % interface between Line 1&2
% row=8:10; start_frame=154; %Line 17
% col=175; 

% FullBondLengths=zeros(layer,numFilaments-1);
% AllBondLengths=zeros(NumTest*layer,numFilaments-1);

% Temp.S1 = load('..\Experiment_Analytical_Abaqus\T240S42\t240s42_1_temps.mat');
% temps = Temp.(strcat('S', int2str(iii)));

temps = load('..\Experiment_Analytical_Abaqus\T240S42\t240s42_1_temps.mat');

% Note: numFilaments-1 = size(row,2)
pixeltemp = zeros(size(row,2), size(temps.(strcat('layer', int2str(layer))),3));
Int_temps_exp = zeros(size(row,2), size(temps.(strcat('layer', int2str(layer))),3));
FrameIdx=zeros(1,numFilaments);


% initialize bond lengths array
Exp_BondLengths=zeros(1,numFilaments-1);
Model_BondLengths=zeros(1,numFilaments-1);

% 1st non-zero temp. data column index for the 1st interface of the experiment
firstEl=1;
% Neck growth calculation using the temperature field from experiment

% % 1st layer from experiment
% jj=1;
    
for rr=1:numFilaments-1 % loop over each interface for layer 1
    
    kk=1;
    % convert Celcius to Kelvin, 174:180 -> find the max temp. near mid
    % point, col=175 is not always the mid point as we go to next lines bc. of
    % angle of the camera
    pixeltemp(rr, :) = max(temps.(strcat('layer', int2str(layer)))(row(rr), 174:180, :))+Kel;
    if rr==1
        [~,sortIdx]=sort(pixeltemp(rr, :));
        FrameIdx(rr)=sortIdx(end);
    else
        [~,sortIdx]=sort(pixeltemp(rr, :));
        while (sortIdx(end-kk)-FrameIdx(rr-1))<2 % make sure that adjacent frames are not taken
            kk=kk+1;
        end
        FrameIdx(rr) = sortIdx(end-kk);
    end
    % interface temperatures at point B
    Int_temps_exp(rr, FrameIdx(rr)-(FrameIdx(1)-1):end-(FrameIdx(1)-1)) = pixeltemp(rr, FrameIdx(rr):end);

    start_frame = 1;
    final_frame = size(Int_temps_exp,2);
    freq = 10.5; % frequency of frames
    exp_time = start_frame/freq:1/freq:final_frame/freq;

    tfinal=exp_time(end);
%     Int_temps_exp( :, ~any(Int_temps_exp,1) ) = [];
    num=size(Int_temps_exp,2)-FrameIdx(rr)+1;
    dt=tfinal/num;
    
    if rr~=1
        % update starting time index (col. no) for that interface, starting from 2nd
        firstEl=firstEl+FrameIdx(rr)-FrameIdx(rr-1);
    end
    
    % interface temperature for experiment
    T_exp = Int_temps_exp(rr,firstEl:end);
    
    % ==================================================================== %
    %% Experimental Bond Length
    % compute Experimental bond length for the 1st layer
    Exp_BondLengths(rr) = Neck_Growth_Model(T_exp,inputs,ProcessParam,...
                    Ref_Temp,tfinal,dt,exp_time,firstEl);
    

end


% 2nd layer from model
jj=2;
    

for rr=1:numFilaments-1 % loop over each interface for layer 1
   
    % ==================================================================== %
    %% Temperature Model
    % nozzle temperature and printer speed for each sample
    TV=[240 42];
    
    % the nozzle temperature is observed to be at least 10 degrees less
    T_N = TV(1)-10; v_p = TV(2);
    
    [Temp_C, time] = FDM_TempModel(Matrix,x,T_N,v_p,ProcessParam);
    
    % convert Celcius to Kelvin
    Temp = Temp_C + Kel;
    
    % Average two adjacent lines' temperatures to obtain interface temperatures
    % for each interface on that layer
%     int_Temp = zeros(size(Temp,1),numFilaments-1);

    ii=rr+1;

    if ii==2
%         kk_int=(jj-1)*(numFilaments-1)+ii;
        kk=(jj-1)*(numFilaments-1)+ii+jj-2;
        T_idx = find(Temp(:,kk)<T_N+Kel, 1, 'first');
        % find the first index that temp value is lower thatn 145 celcius
        T_idx_ = find((Temp(:,kk)+Temp(:,kk-1))/2<Ref_Temp+Kel, 1, 'first');
        int_Temp(1:T_idx_-T_idx+1,rr) = ...
            (Temp(T_idx:T_idx_,kk-1)+...
            Temp(T_idx:T_idx_,kk))/2;
    elseif ii==numFilaments
%         kk_int=(jj-1)*(numFilaments-1)+ii;
        kk=(jj-1)*(numFilaments-1)+ii+jj-1;
        T_idx = find(Temp(:,kk)<T_N+Kel, 1, 'first');
        % find the first index that temp value is lower thatn 145 celcius
        T_idx_ = find((Temp(:,kk)+Temp(:,kk-1))/2<Ref_Temp+Kel, 1, 'first');
        int_Temp(1:T_idx_-T_idx+1,rr) = ...
            (Temp(T_idx:T_idx_,kk-1)+...
            Temp(T_idx:T_idx_,kk))/2;
    else
        kk=(jj-1)*(numFilaments-1)+ii;
        T_idx = find(Temp(:,kk+1)<T_N+Kel, 1, 'first');
        % find the first index that temp value is lower thatn 145 celcius
        T_idx_ = find((Temp(:,kk+1)+Temp(:,kk))/2<Ref_Temp+Kel, 1, 'first');
        int_Temp(1:T_idx_-T_idx+1,rr) = ...
            (Temp(T_idx:T_idx_,kk-1)+...
            Temp(T_idx:T_idx_,kk))/2;
    end

    
    int_Temp_model = int_Temp'; % interface temperature
    
    tfinal=time(T_idx_-T_idx+1);
%     num=size(int_Temp_model,1)-T_idx(1)+1;
    num=size(int_Temp_model,2);
    dt=tfinal/num;
    
    % interface temperature for model, 2nd layer jj=2
    T = int_Temp_model(rr,:);
%     T( :, ~any(T,1) ) = [];



    %% Model Bond Length
    % compute Experimental bond length for the 1st layer
    Model_BondLengths(rr) = Neck_Growth_Model(T,inputs,ProcessParam,...
                    Ref_Temp,tfinal,dt,time,0);
    

end





% FullBondLengths(:,:) = NeckRadius;
% 
% AllBondLengths((iii-1)*layer+1:iii*layer,:) = FullBondLengths;



function[BondLengths] = Neck_Growth_Model(T,inputs,ProcessParam,Ref_temp,tfinal,dt,time,firstEl)

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
    b2 = inputs(4); % model parameter for temp dependent viscosity
    
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
    
    
    




function[temp,abcissa] = FDM_TempModel(matrix,x,T_n,V_p,ProcessParam)

%Process Variables
T_L = T_n; %Extrusion temperature (ºC)
T_E = ProcessParam(1); %Temperature of the envelope (ºC)
v = V_p; % printer speed/Velocity of the extrusion head (m/sec)
%Filament dimensions
Wid = ProcessParam(2); % Layer Thickness (meters)
Ht = ProcessParam(3); % Layer Height (meters)
L = ProcessParam(4); % Length of the filament (meters)

roC = 9.6228*1e5; % Density (kg/m^3)
% h_conv = MaterialParam(2); % Heat transfer coefficient (loss of heat by natural convection)
% lam = MaterialParam(3); % filament and adjacent filament & support
% N_Samples = size(ro,1);

                
%____________________________________ STEP 1 ____________________________________
%Definition of the vector that contains the number of total filaments in each layer
matrix_lin = size(matrix,1);
matrix_col = size(matrix,2);
vector = zeros(matrix_lin,2);
ctr = 0; % # of layers
for i = matrix_lin:-1:1
    ctr = ctr + 1;
    for j = 1:matrix_col
        if matrix(i,j) ~= 0
            vector(ctr,1) = vector(ctr,1) + 1;
        end
    end
end

%Number of layers
m = length(vector(:,1));
%Number of filaments
n = 0;
for j = 1:m
    if m == 1
        n = vector(1,1);
    else
        if vector(j,2) <= 1
            n = n + vector(j,1);
        end
    end
end
%____________________________________ STEP 2 ____________________________________
%Computation variables
dt = .1; %Step time
temp_extra = 10; %Additional time computation after construction of the part
erro = 1e-3; %Convergence error
% FinalTimeFirstFila = .1; % delay time between filamentss
%____________________________________ STEP 3 ____________________________________
%Definition of the size of the variables
% ii = 5; # of contacts with adjacent filaments and support
h = zeros(1,5); lambda = zeros(1,5); a = zeros(n,5); T = zeros (n,5);
vec_b = zeros(n,5); vec_Q = zeros(n,5); b = zeros(1,n); Q = zeros(1,n);
T_begin = zeros(1,n); dif = zeros(1,n); Biot = zeros(1,n); save_T = zeros(1,n);
old_T = zeros(1,n); save_lim = zeros(1,n); viz = zeros(11,n);
%____________________________________ STEP 4 ____________________________________
%Process Variables
for lin = 1:n %Temperature of support (ºC)
    T(lin,5) = T_E;
end

% % For Cylinder
% area = pi * (Wid/2)^2; %Area of the cross section of filament (meters^2) -
% per = pi * Wid; %Perimeter of the cross section of filament (meters)

% For Ellipse
aa=Wid/2; bb=Ht/2;
hh = (aa-bb)^2/(aa+bb)^2;
area = pi * (aa*bb); %Ellipse:Area of the cross section of filament (meters^2)
per = pi * (aa+bb)*(1+3*hh/(10+sqrt(4-3*hh))); %Perimeter of the cross section of filament (meters)

vol = area*L; %Volume of the filament
A_p = per*L; %Superficial area of the filament

% Material Properties
%Thermal conductivity (W/m.K)
% conductivity(1) = 0.1768; % material A
% conductivity(1) = 0.15; % material A
% conductivity(2) = 0.5; % material B
conductivity = 0.15; % material A
% %Density (kg/m^3)
% ro(1) = 1050; % material A
% ro(1) = 1040; % material A
% ro(2) = 1500; % material B
% %Specific heat (J/kg.K)
% C(1) = 2019.7; % material A
% C(1) = 1290; % material A
% C(2) = 2500.7; % material B
%____________________________________ STEP 5 ____________________________________
% Heat transfer coefficient (loss of heat by natural convection)
h_conv = 86;
%Thermal contact conductances between
h(1,1) = 200; % filament and left adjacent filament
h(1,2) = 200; % filament and down adjacent filament
h(1,3) = 200; % filament and right adjacent filament
h(1,4) = 200; % filament and top adjacent filament
h(1,5) = 86; % filament and support
% h(1,5) = h_cond;
%Fraction of perimeter contact between
lambda(1,1) = 0.15; % filament and left adjacent filament
lambda(1,2) = 0.15; % filament and down adjacent filament
lambda(1,3) = 0.15; % filament and right adjacent filament
lambda(1,4) = 0.15; % filament and top adjacent filament
lambda(1,5) = 0.15; % filament and support
% lambda(1,1:5)=lam;
%____________________________________ STEP 6 ____________________________________
%Definition of the parameters influenced by the contacts
for col = 1:5
    for lin = 1:n
        vec_b(lin,col) = h(1,col)*lambda(1,col);
        vec_Q(lin,col) = vec_b(lin,col)*T(lin,col);
    end
end
% %____________________________________ STEP 7 ____________________________________
% %Definition of the parameters influenced by the material properties

scalar = -per/(roC*area);
kt = conductivity;

limite=zeros(n+2,2); % BK

for kk = 1:1:(n+2)
    if mod(kk,2) == 1
        limite(kk,1) = (kk*L-x)/v;
        limite(kk,2) = (kk*L+x)/v;
    else
        limite(kk,1) = limite(kk-1,2);
        limite(kk,2) = ((kk+1)*L-x)/v;
    end
end

% temp=zeros(); % BK
for road = 1:n
    line = 0;
    for i = 0:dt:limite(n,2)
        line = line + 1;
        temp(line,road) = T_L;
    end
end

% abcissa=zeros((limite(end,2)+temp_extra)/dt,1); % BK
for layer = 1:m
    if layer == 1
        for num = 1:vector(layer,1)
            if num == 1
                %____________________________________ STEP 9 ____________________________________
                a(num,5) = 1; %Activation of the contact with support
                %____________________________________ STEP 10 ____________________________________
                %Definition of the variables b and Q defined in equation Eq. 7
                b(num) = h_conv*(1-lambda*a(num,:)') + vec_b(num,:)*a(num,:)';
                Q(num) = (h_conv*(1-lambda*a(num,:)')*T_E + ...
                    vec_Q(num,:)*a(num,:)')/b(num);
                
                %____________________________________ STEP 11 ____________________________________
                p = 0;
                for t = 0:dt:limite(num,1)
                    p = p+1; abcissa(p) = t;
                end
                %____________________________________ STEP 12 ____________________________________
                %Computation of the temperatures of the first filament
                %                 for t = (limite(num,1)+dt):dt:limite_final
                %                 for t = (limite(num,1)+dt):dt:FinalTimeFirstFila % BK
                for t = (limite(num,1)+dt):dt:limite(1,2) % BK
                    p = p+1; abcissa(p) = t;
                    temp(p,num)=(T_L-Q(num))*exp(scalar(num)*b(num)* ...
                        (t-limite(num,1))) + Q(num);
                end
                %Saving the last temperature of the period time of cooling down
                T_begin(num) = temp(p,num);
                %____________________________________ STEP 13 ____________________________________
                %Verification of the value of Biot Number
                Biot(num) = (vol/A_p)*(b(num)/kt(num));
                if Biot(num)>=0.1
                    'WARNING! We cannot use a Lumped System';
                end
                
                %____________________________________ STEP 14 ____________________________________
            else
                %Activation of the contacts
                a(num-1,3) = 1; a(num,1) = 1; a(num,5) = 1;
                %____________________________________ STEP 15 ____________________________________
                %Up-dating of the variable b
                for j = 1:num
                    b(j) = h_conv*(1-lambda*a(j,:)') + vec_b(j,:)*a(j,:)';
                end
                
                %____________________________________ STEP 16 ____________________________________
                if m == 1
                    if num == vector(layer,1)
                        limite_final = limite(num,2) + temp_extra;
                    else
                        limite_final = limite(num,2);
                    end
                else
                    limite_final = limite(num,2);
                end
                
                %                 for t = (FinalTimeFirstFila+dt):dt:limite_final % BK
                for t = (limite(num,1)+dt):dt:limite_final
                    p = p+1; abcissa(p) = t;
                    last = p-1;
                    for j = 1:num
                        save_T(j) = temp(last,j);
                    end
                    %____________________________________ STEP 17 ____________________________________
                    %Iterative process
                    for q = 1:100000
                        %Saving contacts and temperatures of adjacent filaments
                        for j = 1:num
                            if j == 1
                                T(j,3) = save_T(j+1);
                                viz(3,j) = j+1;
                            end
                            if j > 1 && j < num
                                T(j,1) = save_T(j-1);
                                viz(1,j) = j-1;
                                T(j,3) = save_T(j+1);
                                viz(3,j) = j+1;
                            end
                            if j == num
                                T(j,1) = save_T(j-1);
                                viz(1,j) = j-1;
                            end
                            for k = 1:5
                                if T(j,k) ~= 0 && k ~= 5
                                    vec_Q(j,k) = vec_b(j,k)*T(j,k);
                                end
                            end
                            %Up-dating of the variable Q
                            Q(j) = (h_conv*(1-lambda*a(j,:)')*T_E + ...
                                vec_Q(j,:)*a(j,:)')/b(j);
                            old_T(j) = save_T(j);
                        end
                        %Computation of the temperatures
                        if num == 2
                            save_T(num-1) = (T_begin(num-1)-Q(num-1))* ...
                                exp(scalar*b(num-1)*(t-limite(num,1)))+Q(num-1);
                            save_T(num) = (T_L-Q(num))* ...
                                exp(scalar*b(num)*(t-limite(num,1)))+ Q(num);
                            save_lim(1,1) = limite(num,1);
                            save_lim(1,2) = limite(num,1);
                        else
                            for j=1:num-2
                                save_T(j) = (T_begin(j)-Q(j))*exp(scalar*b(j)* ...
                                    (t-limite(j,1)))+Q(j);
                            end
                            save_T(num-1) = (T_begin(num-1)-Q(num-1))* ...
                                exp(scalar*b(num-1)*(t-limite(num,1)))+Q(num-1);
                            save_T(num) = (T_L-Q(num))* ...
                                exp(scalar*b(num)*(t-limite(num,1)))+ Q(num);
                            save_lim(1,num-1) = limite(num,1);
                            save_lim(1,num) = limite(num,1);
                        end
                        for j = 1:num
                            dif(j) = abs(save_T(j)-old_T(j));
                        end
                        try_ = 1;
                        stop = 0;
                        for j = 1:num
                            if dif(try_) < erro
                                try_ = try_+1;
                            end
                            if try_ == num+1
                                stop = 1;
                            end
                        end
                        if stop == 1
                            for j = 1:num
                                temp(p,j) = save_T(j);
                            end
                            break;
                        end
                    end
                end
                T_begin(num) = temp(p,num);
                %End of iterative process
                %____________________________________ STEP 18 ____________________________________
                %Verification of the Biot Number
                for j=1:num
                    Biot(j) = (vol/A_p)*(b(j)/kt);
                    if Biot(j)>=0.1
                        'WARNING! We can not use a Lumped System';
                        j;
                        Biot(j);
                    end
                end
            end
        end
        % =================================================================== %
        % For the remaining layers -> BK
    else
        jj=1;
        for num = sum(vector(1:layer-1,1))+1:sum(vector(1:layer-1,1))+vector(layer,1)
            
            % =================================================================== %
            if mod(layer,2)==0 % for even layers
                
                % =================================================================== %
                % if the filament is the first filament on that layer
                if num == sum(vector(1:layer-1,1))+1
                    %____________________________________ STEP 9 ____________________________________
                    a(num,2) = 1; %Activation of the contact down filament
                    a(num-1,4) = 1; %Activation of the contact top filament
                    
                    % rest of the filaments on that layer
                else % if the filament is not the first filament on that layer
                    jj=jj+1;
                    
                    %____________________________________ STEP 14 ____________________________________
                    %Activation of the contacts (horizontal)
                    a(num-1,1) = 1; a(num,3) = 1;
                    %Activation of the contacts (vertical)
                    a(num,2) = 1; a(num-(2*jj-1),4) = 1;
                end
                
                %____________________________________ STEP 15 ____________________________________
                %Up-dating of the variable b
                for j = 1:num
                    b(j) = h_conv*(1-lambda*a(j,:)') + vec_b(j,:)*a(j,:)';
                end
                
                %____________________________________ STEP 16 ____________________________________
                if m == layer
                    if num == sum(vector(:,1))
                        limite_final = limite(num,2) + temp_extra;
                    else
                        limite_final = limite(num,2);
                    end
                else
                    limite_final = limite(num,2);
                end
                
                % check is array element is empty or not
                chk = size(save_T,2);
                
                for t = (limite(num,1)+dt):dt:limite_final
                    p = p+1; abcissa(p) = t;
                    last = p-1;
                    for j = 1:num
                        save_T(j) = temp(last,j);
                    end
                    
                    %____________________________________ STEP 17 ____________________________________
                    %Iterative process
                    for q = 1:10000
                        %Saving contacts and temperatures of adjacent filaments
                        for j = 1:num
                            if j == 1
                                T(j,3) = save_T(j+1);
                                viz(3,j) = j+1;
                                if chk >= j+2*(sum(vector(1,1))-j)+1
                                    T(j,4) = save_T(j+2*(sum(vector(1,1))-j));
                                    viz(4,j) = j+2*(sum(vector(1,1))-j);
                                end
                            end
                            if j > 1 && j < sum(vector(1,1))
                                T(j,1) = save_T(j-1);
                                viz(1,j) = j-1;
                                T(j,3) = save_T(j+1);
                                viz(3,j) = j+1;
                                if chk >= j+2*(sum(vector(1,1))-j)+1
                                    T(j,4) = save_T(j+2*(sum(vector(1,1))-j));
                                    viz(4,j) = j+2*(sum(vector(1,1))-j);
                                end
                            end
                            if j == sum(vector(1,1))
                                T(j,1) = save_T(j-1);
                                viz(1,j) = j-1;
                                if chk >= j+1
                                    T(j,4) = save_T(j+1);
                                    viz(4,j) = j+1;
                                end
                            end
                            
                            % loop over all even layers except first one
                            for kk=2:2:m
                                % even layers' first filament
                                if j == sum(vector(1:kk-1,1))+1
                                    T(j,1) = save_T(j+1);
                                    viz(1,j) = j+1;
                                    T(j,2) = save_T(j-1);
                                    viz(2,j) = j-1;
                                end
                                if j > sum(vector(1:kk-1,1))+1 && j < sum(vector(1:kk-1,1))+vector(kk,1)
                                    T(j,3) = save_T(j-1);
                                    viz(3,j) = j-1;
                                    T(j,1) = save_T(j+1);
                                    viz(1,j) = j+1;
                                    if chk >= j-(2*(j-(sum(vector(1:kk-1,1))+1))+1)
                                        T(j,2) = save_T(j-(2*(j-(sum(vector(1:kk-1,1))+1))+1));
                                        viz(2,j) = j-(2*(j-(sum(vector(1:kk-1,1))+1))+1);
                                    end
                                end
                                if j == sum(vector(1:kk-1,1))+vector(kk,1)
                                    T(j,3) = save_T(j-1);
                                    viz(3,j) = j-1;
                                    if chk >= j-(2*(j-(sum(vector(1:kk-1,1))+1))+1)
                                        T(j,2) = save_T(j-(2*(j-(sum(vector(1:kk-1,1))+1))+1));
                                        viz(2,j) = j-(2*(j-(sum(vector(1:kk-1,1))+1))+1);
                                    end
                                end
                            end
                            
                            % loop over all odd layers except first one
                            for kk=3:2:m
                                % odd layers' first filament
                                if j == sum(vector(1:kk-1,1))+1
                                    T(j,3) = save_T(j+1);
                                    viz(3,j) = j+1;
                                    T(j,2) = save_T(j-1);
                                    viz(2,j) = j-1;
                                    T(j,4) = save_T(j+(2*(j-(sum(vector(1:kk-1,1))+1))));
                                    viz(4,j) = j+(2*(j-(sum(vector(1:kk-1,1))+1)));
                                end
                                % odd layers' inner filaments
                                if j > sum(vector(1:kk-1,1))+1 && j < sum(vector(1:kk-1,1))+vector(kk,1)
                                    T(j,1) = save_T(j-1);
                                    viz(1,j) = j-1;
                                    T(j,3) = save_T(j+1);
                                    viz(3,j) = j+1;
                                    T(j,2) = save_T(j-(2*(j-(sum(vector(1:kk-1,1))+1))+1));
                                    viz(2,j) = j-(2*(j-(sum(vector(1:kk-1,1))+1))+1);
                                    T(j,4) = save_T(j+(2*(j-(sum(vector(1:kk,1))+2))));
                                    viz(4,j) = j+(2*(j-(sum(vector(1:kk-1,1))+2)));
                                end
                                % odd layers' last filament
                                if j == sum(vector(1:kk-1,1))+vector(kk,1)
                                    T(j,1) = save_T(j-1);
                                    viz(1,j) = j-1;
                                    T(j,2) = save_T(j-(2*(j-(sum(vector(1:kk-1,1))+1))+1));
                                    viz(2,j) = j-(2*(j-(sum(vector(1:kk-1,1))+1))+1);
                                    T(j,4) = save_T(j+1);
                                    viz(4,j) = j+1;
                                end
                            end
                            
                            for k = 1:5
                                if T(j,k) ~= 0 && k ~= 5
                                    vec_Q(j,k) = vec_b(j,k)*T(j,k);
                                end
                            end
                            %Up-dating of the variable Q
                            Q(j) = (h_conv*(1-lambda*a(j,:)')*T_E + ...
                                vec_Q(j,:)*a(j,:)')/b(j);
                            old_T(j) = save_T(j);
                        end
                        
                        %Computation of the temperatures
                        for j=1:num-2
                            save_T(j) = (T_begin(j)-Q(j))*exp(scalar*b(j)* ...
                                (t-limite(j,1)))+Q(j);
                        end
                        save_T(num-1) = (T_begin(num-1)-Q(num-1))*exp(scalar*b(num-1)* ...
                            (t-limite(num,1)))+Q(num-1);
                        save_T(num) = (T_L-Q(num))*exp(scalar*b(num)*(t-limite(num,1)))+Q(num);
                        save_lim(1,num-1) = limite(num,1);
                        save_lim(1,num) = limite(num,1);
                        
                        for j = 1:num
                            dif(j) = abs(save_T(j)-old_T(j));
                        end
                        try_ = 1;
                        stop = 0;
                        for j = 1:num
                            if dif(try_) < erro
                                try_ = try_+1;
                            end
                            if try_ == num+1
                                stop = 1;
                            end
                        end
                        if stop == 1
                            for j = 1:num
                                temp(p,j) = save_T(j);
                            end
                            break;
                        end
                    end
                end
                T_begin(num) = temp(p,num);
                %End of iterative process
                %____________________________________ STEP 18 ____________________________________
                %Verification of the Biot Number
                for j=1:num
                    Biot(j) = (vol/A_p)*(b(j)/kt);
                    if Biot(j)>=0.1
                        'WARNING! We can not use a Lumped System';
                        j;
                        Biot(j);
                    end
                end
                
                % =================================================================== %
            else % odd layers
                % =================================================================== %
                % if the filament is the first filament on that layer
                if num == sum(vector(1:layer-1,1))+1
                    %____________________________________ STEP 9 ____________________________________
                    a(num,2) = 1; %Activation of the contact down filament
                    a(num-1,4) = 1; %Activation of the contact top filament
                    
                    % rest of the filaments on that layer
                else % if the filament is not the first filament on that layer
                    jj=jj+1;
                    
                    %____________________________________ STEP 14 ____________________________________
                    %Activation of the contacts (horizontal)
                    a(num-1,3) = 1; a(num,1) = 1;
                    %Activation of the contacts (vertical)
                    a(num,2) = 1; a(num-(2*jj-1),4) = 1;
                end
                
                %____________________________________ STEP 15 ____________________________________
                %Up-dating of the variable b
                for j = 1:num
                    b(j) = h_conv*(1-lambda*a(j,:)') + vec_b(j,:)*a(j,:)';
                end
                
                %____________________________________ STEP 16 ____________________________________
                if m == layer
                    if num == sum(vector(:,1))
                        limite_final = limite(num,2) + temp_extra;
                    else
                        limite_final = limite(num,2);
                    end
                else
                    limite_final = limite(num,2);
                end
                
                % check is array element is empty or not
                chk = size(save_T,2);
                
                for t = (limite(num,1)+dt):dt:limite_final
                    p = p+1; abcissa(p) = t;
                    last = p-1;
                    for j = 1:num
                        save_T(j) = temp(last,j);
                    end
                    
                    %____________________________________ STEP 17 ____________________________________
                    %Iterative process
                    for q = 1:100000
                        %Saving contacts and temperatures of adjacent filaments
                        for j = 1:num
                            %Saving contacts and temperatures of adjacent filaments
                            if j == 1
                                T(j,3) = save_T(j+1);
                                viz(3,j) = j+1;
                                T(j,4) = save_T(j+2*(sum(vector(1,1))-j));
                                viz(4,j) = j+2*(sum(vector(1,1))-j);
                            end
                            if j > 1 && j < sum(vector(1,1))
                                T(j,1) = save_T(j-1);
                                viz(1,j) = j-1;
                                T(j,3) = save_T(j+1);
                                viz(3,j) = j+1;
                                T(j,4) = save_T(j+2*(sum(vector(1,1))-j));
                                viz(4,j) = j+2*(sum(vector(1,1))-j);
                            end
                            if j == sum(vector(1,1))
                                T(j,1) = save_T(j-1);
                                viz(1,j) = j-1;
                                T(j,4) = save_T(j+1);
                                viz(4,j) = j+1;
                            end
                            
                            % loop over all even layers except first one
                            for kk=2:2:m
                                % even layers' first filament
                                if j == sum(vector(1:kk-1,1))+1
                                    T(j,1) = save_T(j+1);
                                    viz(1,j) = j+1;
                                    T(j,2) = save_T(j-1);
                                    viz(2,j) = j-1;
                                    T(j,4) = save_T(j+(2*((sum(vector(1:kk-1,1))-j))));
                                    viz(4,j) = j+(2*((sum(vector(1:kk-1,1))-j)));
                                end
                                if j > sum(vector(1:kk-1,1))+1 && j < sum(vector(1:kk-1,1))+vector(kk,1)
                                    T(j,1) = save_T(j+1);
                                    viz(1,j) = j+1;
                                    T(j,2) = save_T(j-(2*(j-(sum(vector(1:kk-1,1))+1))+1));
                                    viz(2,j) = j-(2*(j-(sum(vector(1:kk-1,1))+1))+1);
                                    T(j,3) = save_T(j-1);
                                    viz(3,j) = j-1;
                                    %                                     j+(2*((sum(vector(1:kk,1))-j))+1)
                                    T(j,4) = save_T(j+(2*((sum(vector(1:kk-1,1))-j)))+1);
                                    viz(4,j) = j+(2*((sum(vector(1:kk-1,1))-j)));
                                end
                                if j == sum(vector(1:kk-1,1))+vector(kk,1)
                                    T(j,3) = save_T(j-1);
                                    viz(3,j) = j-1;
                                    T(j,4) = save_T(j+1);
                                    viz(4,j) = j+1;
                                    T(j,2) = save_T(j-(2*(j-(sum(vector(1:kk-1,1))+1))+1));
                                    viz(2,j) = j-(2*(j-(sum(vector(1:kk-1,1))+1))+1);
                                end
                            end
                            
                            % loop over all odd layers except first one
                            for kk=3:2:m
                                % odd layers' first filament
                                if j == sum(vector(1:kk-1,1))+1
                                    T(j,3) = save_T(j+1);
                                    viz(3,j) = j+1;
                                    T(j,2) = save_T(j-1);
                                    viz(2,j) = j-1;
                                end
                                % odd layers' inner filaments
                                if j > sum(vector(1:kk-1,1))+1 && j < sum(vector(1:kk-1,1))+vector(kk,1)
                                    T(j,1) = save_T(j-1);
                                    viz(1,j) = j-1;
                                    T(j,3) = save_T(j+1);
                                    viz(3,j) = j+1;
                                    if chk >= j-(2*(j-(sum(vector(1:kk-1,1))+1))+1)
                                        T(j,2) = save_T(j-(2*(j-(sum(vector(1:kk-1,1))+1))+1));
                                        viz(2,j) = j-(2*(j-(sum(vector(1:kk-1,1))+1))+1);
                                    end
                                end
                                % odd layers' last filament
                                if j == sum(vector(1:kk-1,1))+vector(kk,1)
                                    T(j,1) = save_T(j-1);
                                    viz(1,j) = j-1;
                                    if chk >= j-(2*(j-(sum(vector(1:kk-1,1))+1))+1)
                                        T(j,2) = save_T(j-(2*(j-(sum(vector(1:kk-1,1))+1))+1));
                                        viz(2,j) = j-(2*(j-(sum(vector(1:kk-1,1))+1))+1);
                                    end
                                end
                            end
                            
                            for k = 1:5
                                if T(j,k) ~= 0 && k ~= 5
                                    vec_Q(j,k) = vec_b(j,k)*T(j,k);
                                end
                            end
                            %Up-dating of the variable Q
                            Q(j) = (h_conv*(1-lambda*a(j,:)')*T_E + ...
                                vec_Q(j,:)*a(j,:)')/b(j);
                            old_T(j) = save_T(j);
                        end
                        %Computation of the temperatures
                        for j=1:num-2
                            save_T(j) = (T_begin(j)-Q(j))*exp(scalar*b(j)* ...
                                (t-limite(j,1)))+Q(j);
                        end
                        save_T(num-1) = (T_begin(num-1)-Q(num-1))*exp(scalar*b(num-1)* ...
                            (t-limite(num,1)))+Q(num-1);
                        save_T(num) = (T_L-Q(num))*exp(scalar*b(num)*(t-limite(num,1)))+Q(num);
                        save_lim(1,num-1) = limite(num,1);
                        save_lim(1,num) = limite(num,1);
                        
                        for j = 1:num
                            dif(j) = abs(save_T(j)-old_T(j));
                        end
                        try_ = 1;
                        stop = 0;
                        for j = 1:num
                            if dif(try_) < erro
                                try_ = try_+1;
                            end
                            if try_ == num+1
                                stop = 1;
                            end
                        end
                        if stop == 1
                            for j = 1:num
                                temp(p,j) = save_T(j);
                            end
                            break;
                        end
                    end
                end
                
                T_begin(num) = temp(p,num);
                %End of iterative process
                %____________________________________ STEP 18 ____________________________________
                %Verification of the Biot Number
                for j=1:num
                    Biot(j) = (vol/A_p)*(b(j)/kt);
                    if Biot(j)>=0.1
                        'WARNING! We can not use a Lumped System';
                        j;
                        Biot(j);
                    end
                end
                
            end
        end
    end
end


end






% end


%     
%     % Load Gaussian Process model for K1 and Work done
%     GPStruc_K1 = load('06 GP Model - K1 032 5-90\Fitting\gprMdl_K1');
%     GPStruc_W = load('06 GP Model - W 032 5-90\Fitting\gprMdl_W');
%     GP = {GPStruc_K1 GPStruc_W};
% 
%     %% =================   surrogateopt Inputs  ===========================
%     a_crit = 19.337218673249; % mm
%     
%     % randomly choose Minimum Work done obtained from global opt.
%     W_Global = 4.232302423175552e+04;
%     t_global = 4219; % [s]
%     
%     % mission 1 work done assuming 4 missions, dividing into 4 intervals
% %     W_total_1 = W_Global*0.30; % 1st mission 30% of total work done
% %     t_total_1 = t_global*0.30; % 1st mission 30% of total # of cycles
% %     W_total_1 = W_Global*0.20; % 2nd mission 20% of total work done
% %     t_total_1 = t_global*0.20; % 2nd mission 20% of total # of cycles
% %     W_total_1 = W_Global*0.40; % 3rd mission 40% of total work done
% %     t_total_1 = t_global*0.40; % 3rd mission 40% of total # of cycles
%     W_total_1 = W_Global*0.10; % 3rd mission 40% of total work done
%     t_total_1 = t_global*0.10; % 3rd mission 40% of total # of cycles
%     
%     % lower and upper bound constraints
%     F1_lb=3000; F1_ub=4000;
%     F2_lb=5000; F2_ub=5400;
% 	F3_lb=4000; F3_ub=5000;
% 
% %     lb = [F1_lb F2_lb F3_lb 295 395 295];
% %     ub = [F1_ub F2_ub F3_ub 305 405 305];
%     lb = [F1_lb F2_lb F3_lb 0.2*t_total_1*0.8 0.6*t_total_1*0.8 0.2*t_total_1*0.8];
%     ub = [F1_ub F2_ub F3_ub 0.2*t_total_1*1.2 0.6*t_total_1*1.2 0.2*t_total_1*1.2];
% 
%     % 1st mission appx 5000 cycles
%     
%     % fixed frequencies of each block loading
%     w1 = 5; w2 = 5; w3 = 5;
%     w = [w1 w2 w3];
% 	
%     % Initial value for design variables; amplitude, time
%     % [A1 A2 A3 t1 t2 t3]
%     x_initial = [3500 5200 4500 0.2*t_total_1 0.6*t_total_1 0.2*t_total_1];
% 
%     %% MC Samples
%     % % Mean of C, material parameter
%     mu_C = 4.5628e-9; cov_C = 0.10;
%     
%     % icdf of the parameters
%     U = lhsdesign(n_samples, 1);
%     C_Samples = icdf('normal', U(:, 1), mu_C, mu_C*cov_C)
% 
% 
% %     W_total_2 = W_Global*0.30; % 1st mission 30% of total work done
% %     W_total_3 = W_Global*0.40; % 1st mission 40% of total work done
% %     W_total_4 = W_Global*0.10; % 1st mission 10% of total work done
%     %     W_total_1 = rand*(W_Global/4 - 0) + 0
% %     W_total_2 = rand*(W_Global/2 - W_Global/4) + W_Global/4;
%     
%     objectivefun = @(x)objfun(x,w,a0,m,n_samples,C_Samples,GP,W_total_1,t_total_1,a_crit);
%                           
%      options = optimoptions('surrogateopt','PlotFcn','surrogateoptplot',...
%         'InitialPoints',x_initial,'UseParallel',true);%,'MaxFunctionEvaluations',30);
%     %,'ObjectiveLimit',1/n_samples-1e-13);
% 
%     % Call Optimization
%     [Fsol,Fval,ex,out,trials] = surrogateopt(objectivefun,lb,ub,options);
%     
%     %% %%%%%%%%%%%%%%%%%%     Objective function    %%%%%%%%%%%%%%%%%%%%%%%
%     function [Pf] = objfun(x,w,a0,m,n_samples,C_Samples,GP,W_total_1,t_total_1,a_crit) 
%         
%         % Objective function that needs to be minimized: P[a(Nf) > a_crit]
%         % maximum loadings
%         Fmax_1=x(1); Fmax_2=x(2); Fmax_3=x(3);
%         % min. loadings => R = 0.5
%         Fmin_1=0.5*Fmax_1; Fmin_2=0.5*Fmax_2; Fmin_3=0.5*Fmax_3;
%         
%         % final number of cycles for individual blocks: (freq x time)
%         nf_1 = w(1)*x(4); nf_2 = w(2)*x(5); nf_3 = w(3)*x(6);
%         
% 		% collect into arrays
% 		Fmax = [Fmax_1 Fmax_2 Fmax_3]; Fmin = [Fmin_1 Fmin_2 Fmin_3];
%         N_f = [nf_1 nf_2 nf_3];
% 		
%         % Crack propagation
%         % NumTimes: # of times the predicted final crack size is greater
%         % than the critical crack size
%         [W, a_f, NumTimes] = Forman_K_W_gp(GP,Fmax,Fmin,a0,...
%                                         m,n_samples,N_f,C_Samples,a_crit);
%          
%         % Prob. of failure
%         Pf = (NumTimes)/n_samples;
%         
%         % total duration of 3 blocks
%         t_total = x(4) + x(5) + x(6);
% 
% %         Pf = mean(a_f);
% 
%         % penalty function: total work done in that mission is less than
%         % the minimum required work done that is randomly chosen from the
%         % total work done during the service life, which is calculated
%         % using the global optimization
%         if ( sum(W) <  W_total_1)
%             Pf = 2*Pf;
%         end
%         % penalty function
%         if (t_total < t_total_1) % 20% of the total time 1st mission
%            Pf = 2*Pf;
%         end
% 
%     end
% 
% 
% %% %%%%%%%%%%%%%%%%%%%%        Display & Print     %%%%%%%%%%%%%%%%%%%%%%%%
% % tmp = load('Work.mat','a_final');
% % aF = tmp.('a_final');
% 
% % tmp = load('Work.mat','af');
% % aF = tmp.('af');
% 
% % assignin('base','af',af);
% % assignin('base','a_f',a_f)
% 
% % % Display and Plotting
% fprintf('Amplitude of 1st block: %1.1f [lbs]\n',Fsol(1));
% fprintf('Amplitude of 2nd block: %1.1f [lbs]\n',Fsol(2));
% fprintf('Amplitude of 3rd block: %1.1f [lbs]\n',Fsol(3));
% % fprintf('Frequency of first block: %1.3f [s^-1]\n',Fsol(3));
% % fprintf('Frequency of second block: %1.3f [s^-1]\n',Fsol(4));
% fprintf('Duration of 1st block loading: %.1f [s]\n',Fsol(4));
% fprintf('Duration of 2nd block loading: %.1f [s]\n',Fsol(5));
% fprintf('Duration of 3rd block loading: %.1f [s]\n',Fsol(6));
% fprintf('# of cycles of the 1st block loading: %.1f\n',round(Fsol(4)*w(1)));
% fprintf('# of cycles of the 2nd block loading: %.1f\n',round(Fsol(5)*w(2)));
% fprintf('# of cycles of the 3rd block loading: %.1f\n',round(Fsol(6)*w(3)));
% fprintf('Total # of cycles: %.1f [s]\n',round(w(1)*Fsol(4)+w(2)*Fsol(5)+w(3)*Fsol(6)));
% % fprintf('Probability of failure: %f\n',floor(Fval/a_crit)/n_samples);
% fprintf('Objectove function: %f\n',Fval);
% % if F>=1
% %     fprintf('%d constraint(s) is/(are) not satisfied \n',floor(Fval));
% % else
% %     fprintf('Probability of failure: %f\n',(Fval-floor(Fval)));
% % end
% 
% 
% % Plot block loads    
% % Define some parameters that define the triangle wave.
% elementsPerHalfPeriod1 = Fsol(4)/2; % Number of elements in each rising or falling section.
% elementsPerHalfPeriod2 = Fsol(5)/2; % Number of elements in each rising or falling section.
% elementsPerHalfPeriod3 = Fsol(6)/2; % Number of elements in each rising or falling section.
% % Peak-to-peak amplitude.
% fmax1 = Fsol(1); fmax2 = Fsol(2); fmax3 = Fsol(3); % Peak-to-peak amplitude
% verticalOffset = -2; % Also acts as a phase shift.
% numberOfPeriods1 = w(1); % How many replicates of the triangle you want.
% numberOfPeriods2 = w(2); % Basically the frequencies
% numberOfPeriods3 = w(3);
% % Construct one cycle, up and down.
% risingSignal1 = linspace(0.5*fmax1, fmax1, elementsPerHalfPeriod1);
% fallingSignal1 = linspace(fmax1, 0.5*fmax1, elementsPerHalfPeriod1);
% risingSignal2 = linspace(0.5*fmax2, fmax2, elementsPerHalfPeriod2);
% fallingSignal2 = linspace(fmax2, 0.5*fmax2, elementsPerHalfPeriod2);
% risingSignal3 = linspace(0.5*fmax3, fmax3, elementsPerHalfPeriod3);
% fallingSignal3 = linspace(fmax3, 0.5*fmax3, elementsPerHalfPeriod3);
% % Combine rising and falling sections into one single triangle.
% firstCycle = [risingSignal1, fallingSignal1(2:end)] + verticalOffset; 
% secondCycle = [risingSignal2, fallingSignal2(2:end)] + verticalOffset;
% thirdCycle = [risingSignal3, fallingSignal3(2:end)] + verticalOffset; 
% % Now replicate this cycle several (numberOfPeriods) times.
% waveform1 = repmat(firstCycle, [1 round(numberOfPeriods1)]);
% x1 = 0 : length(waveform1)-1;
% waveform2 = repmat(secondCycle, [1 round(numberOfPeriods2)]);
% x2 = x1(end) : x1(end)+length(waveform2)-1;
% waveform3 = repmat(thirdCycle, [1 round(numberOfPeriods3)]);
% x3 = x2(end) : x2(end)+length(waveform3)-1;
% % Now plot the triangle wave.
% figure();
% n0=27500; % initial # of cycles
% plot(n0+x1,waveform1,'b-',n0+x2,waveform2,'r-',n0+x3,waveform3,'g-','LineWidth',2);
% grid on;
% title('Load Spectrum', 'FontName', 'Times New Roman', ...
%                     'FontSize',18,'Color','k', 'Interpreter', 'LaTeX')
% xlabel('Cycles', 'FontName', 'Times New Roman', ...
%                     'FontSize',18,'Color','k', 'Interpreter', 'LaTeX')
% ylabel('Load [\textit{lbf}]', 'FontName', 'Times New Roman', ...
%                     'FontSize',18,'Color','k', 'Interpreter', 'LaTeX')
% 
% 
%     
%     
%     
% % %% Plot PDF of crack size
% % figure();
% % a_F = a_F * 1e3; % in [mm]
% % a_crit = a_crit*1e3; % in [mm]
% % [fm,xm] = ksdensity(a_F(:));
% % [~, index] = min(abs(xm-a_crit));
% % plot(xm,fm,'.-r')
% % hold on
% % line([a_crit a_crit], [min(fm) max(fm)*(1+0.1) ]);
% % text(a_crit+a_crit*.00003,max(fm),'\rightarrow a_{crit}','FontSize',16);
% % hold on
% % p = xm(index:end);
% % q = fm(index:end);
% % H = area(p,q);
% % set(H(1),'FaceColor',[1 0.5 0])
% % xlabel('Crack size [\textit{mm}]', 'FontName', 'Times New Roman', ...
% %                     'FontSize',18,'Color','k', 'Interpreter', 'LaTeX');
% % ylabel('PDF', 'FontName', 'Times New Roman', ...
% %                     'FontSize',18,'Color','k', 'Interpreter', 'LaTeX');
% % xlim([a0*1e3 a_crit+2e-1])
% 
%     
% 
% end