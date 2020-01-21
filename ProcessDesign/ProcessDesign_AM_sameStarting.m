function [Fsol,Fval,trials] = ProcessDesign_AM_sameStarting(Matrix,x_cut,layerID,ProcessParam,...
                                            inputs_Model,Ref_Temp,TV,...
                                            GPR_model,Coordinates,n_samples)

% convert Celcius to Kelvin
Kel = 273.15; % 0 Celcius in K

% nozzle temperature - 10, b\c the printer does not print at prescribed
% temperature.
T_n = TV(1)-10;
v_p = TV(2);

numFilaments = ProcessParam(4);

Model_BondLengths=zeros(1,numFilaments-1);
% Model_BondLengths_layer1=zeros(1,numFilaments-1);

% lower and upper bound constraints for increment temperature
% boundaries are 10 less than actual b/c printer loses 10 degrees when it
% prints the filament
Tinc_lb=220-10; Tinc_ub=260-10;
vp_lb=0.015; vp_ub=0.043;

lb = [Tinc_lb,vp_lb];
ub = [Tinc_ub,vp_ub];


% Input for the free deposition model (1D heat transfer model)
% alpha = density*Specific heat, it is calibrated together with b2 assuming
% std for measurement error 0.05 mm.

%% MC Samples
% % Mean of C, material parameter
mu_b2 = 0.00379042459;  std_b2 = 0.0037430434;
mu_alpha = 1.1971103278688e6; std_alpha = 3.15969052094e5;

% icdf of the parameters
U = lhsdesign(n_samples, 2);
b2_Samples = icdf('normal', U(:, 1), mu_b2, std_b2);
alpha_Samples = icdf('normal', U(:, 2), mu_alpha, std_alpha);


% Initial value for design variables; increment temperature
x_initial = [230,0.035];
objectivefun = @(x)objfun(x,GPR_model,Coordinates,n_samples,b2_Samples,...
                                            alpha_Samples,T_n,v_p);

options = optimoptions('surrogateopt','PlotFcn','surrogateoptplot',...
    'InitialPoints',x_initial,'UseParallel',false,'MaxFunctionEvaluations',100);
%,'ObjectiveLimit',1/n_samples-1e-13);

% Call Optimization
[Fsol,Fval,ex,out,trials] = surrogateopt(objectivefun,lb,ub,options);


%% %%%%%%%%%%%%%%%%%%     Objective function    %%%%%%%%%%%%%%%%%%%%%%%
    function [Objective] = objfun(x,GPR_model,Coordinates,n_samples,...
                                    b2_Samples,alpha_Samples,T_N,V_p)
        
        % Objective function that needs to be minimized:
%         T_opt=x(1); v_opt=x(2);
        
        % 2nd layer from model
        jj=layerID;

        Mu = zeros(1,n_samples);
        Sigma = zeros(1,n_samples);
        
        % MC Samples loop
        for nn = 1:n_samples
            % loop over each interface for layerID
            for rr=1:numFilaments-1

                % ==================================================================== %
                %% Temperature Model

                [Temp_C, time] = FDM_TempModel(Matrix,x_cut,T_N,V_p,...
                                ProcessParam,alpha_Samples(nn),x);

                % convert Celcius to Kelvin
                Temp = Temp_C + Kel;

                % Average two adjacent lines' temperatures to obtain interface temperatures
                % for each interface on that layer
                %     int_Temp = zeros(size(Temp,1),numFilaments-1);

                ii=rr+1;

                if jj~=1 && ii==2
                    %         kk_int=(jj-1)*(numFilaments-1)+ii;
                    kk=(jj-1)*(numFilaments-1)+ii+jj-2;
                    T_idx = find(Temp(:,kk)<x(1)+Kel, 1, 'first');
                    % find the first index that temp value is lower thatn 145 celcius
                    T_idx_ = find((Temp(:,kk)+Temp(:,kk-1))/2<Ref_Temp+Kel, 1, 'first');
                    int_Temp(1:T_idx_-T_idx+1,rr) = ...
                        (Temp(T_idx:T_idx_,kk-1)+...
                        Temp(T_idx:T_idx_,kk))/2;
                elseif jj~=1 && ii==numFilaments
                    %         kk_int=(jj-1)*(numFilaments-1)+ii;
                    kk=(jj-1)*(numFilaments-1)+ii+jj-1;
                    T_idx = find(Temp(:,kk)<x(1)+Kel, 1, 'first');
                    % find the first index that temp value is lower thatn 145 celcius
                    T_idx_ = find((Temp(:,kk)+Temp(:,kk-1))/2<Ref_Temp+Kel, 1, 'first');
                    int_Temp(1:T_idx_-T_idx+1,rr) = ...
                        (Temp(T_idx:T_idx_,kk-1)+...
                        Temp(T_idx:T_idx_,kk))/2;
                else
                    kk=(jj-1)*(numFilaments-1)+ii+(jj-2);
                    T_idx = find(Temp(:,kk+1)<x(1)+Kel, 1, 'first');
                    % find the first index that temp value is lower thatn 145 celcius
                    T_idx_ = find((Temp(:,kk+1)+Temp(:,kk))/2<Ref_Temp+Kel, 1, 'first');
                    int_Temp(1:T_idx_-T_idx+1,rr) = ...
                        (Temp(T_idx:T_idx_,kk)+...
                        Temp(T_idx:T_idx_,kk+1))/2;
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
                Model_BondLengths(rr) = Neck_Growth_Model(T,inputs_Model,ProcessParam,...
                    Ref_Temp,tfinal,dt,time,0,b2_Samples(nn));

                % gp model input requirement:
                %               T_noz: Nozzle temperature in C,
                %               v_print: printer speed in mm/sec
                %               x and y coordinates
                T_noz = x(1)+10; % add 10 b\c GP model takes input as the 
                % input to the printer, it does not care about the IR
                % camera reading
                v_print = x(2);
                
                
                Gp_input = [T_noz*ones(1*(numFilaments-1),1), ...
                    v_print*ones(1*(numFilaments-1),1), ...
                    Coordinates];
                %__________________________________________________________________________
                
                % predict BL model discrepency using GP inputs for model Temp. data
                [delta_BL_model] = predict(GPR_model, Gp_input);


                % corrected bond lenghs computed using 1D model Temperature data
%                 Corr_Model_BondLengths = Model_BondLengths + ...
%                     delta_BL_model((jj-1)*(numFilaments-2)+jj:jj*(numFilaments-1))';
                Corr_Model_BondLengths = Model_BondLengths + ...
                    delta_BL_model';

            end

            % predicted mean&std of BL for each layer
            Mu(:,nn) = mean(Corr_Model_BondLengths);
            Sigma(:,nn) = std(Corr_Model_BondLengths);
        end
%         Mu = mean(meanModel_2);
%         Std = std(meanModel_2);


        % max. mean of Mu=mean_BL & min. std of Mu=mean_BL
        F1 = -mean(Mu)+std(Mu)

        % min. mean of Sigma=sigma_BL & min. std of Sigma=sigma_BL
        F2 = mean(Sigma)+std(Sigma)


        Objective = F1 / F2;
        
%         % maximize Beta=mu-sigma (better BL with less std/variation)
%         % minimize -mu+sigma
%         % Define Beta = -mu_BL + sigma_BL
%         Beta = -Mu + Sigma;
%         
%         mean(Beta)
%         std(Beta)
%         % minimize Objective
%         weight = 0.5;
%         Objective = weight*mean(Beta) + (1-weight)*std(Beta);
    end

%     % penalty function: total work done in that mission is less than
%     % the minimum required work done that is randomly chosen from the
%     % total work done during the service life, which is calculated
%     % using the global optimization
%     if ( sum(W) <  W_total_1)
%         Pf = 2*Pf;
%     end
%     % penalty function
%     if (t_total < t_total_1) % 20% of the total time 1st mission
%         Pf = 2*Pf;
%     end
end





function[temp,abcissa] = FDM_TempModel(matrix,x,T_n,V_p,ProcessParam,alpha,...
                                OptimVar)

%Process Variables
T_L = T_n; % Extrusion temperature (ºC)
T_E = ProcessParam(5); %Temperature of the envelope (ºC)
v = V_p; % printer speed/Velocity of the extrusion head (m/sec)

% printer speed for the current layer that is optimized process design
T_design = OptimVar(1);
v_design = OptimVar(2);


%Filament dimensions
Wid = ProcessParam(1); % Layer Thickness (meters)
Ht = ProcessParam(2); % Layer Height (meters)
L = ProcessParam(6); % Length of the filament (meters)

roC = alpha; % Density (kg/m^3)

                
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
T_begin = zeros(1,n); dif = zeros(1,n); save_T = zeros(1,n);
% Biot = zeros(1,n); 
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

% vol = area*L; %Volume of the filament
% A_p = per*L; %Superficial area of the filament

% Material Properties
%Thermal conductivity (W/m.K)
% conductivity(1) = 0.1768; % material A
% conductivity(1) = 0.15; % material A
% conductivity(2) = 0.5; % material B
% conductivity = 0.15; % material A
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
% kt = conductivity;

limite=zeros(n+2,2); % BK
kk=0;
for ii = 1:matrix_lin
    kkk=0;
    for jj=1:matrix_col
        if ii==1 && ii~=matrix_lin % 1st layer
            kkk=kkk+1;
            if mod(kkk,2) == 1
                limite(kkk,1) = (kkk*L-x)/v(ii);
                limite(kkk,2) = (kkk*L+x)/v(ii);
            else
                limite(kkk,1) = limite(kkk-1,2);
                limite(kkk,2) = ((kkk+1)*L-x)/v(ii);
            end
        elseif ii==matrix_lin && ii~=1% top/current layer & not 1st layer
            kk=kk+1;
            % odd filaments for even layers, even filaments for odd layers
            if mod(kk,2) == 1
                limite((ii-1)*matrix_col+kk,1) = limite((ii-1)*matrix_col+kk-1,2);
                limite((ii-1)*matrix_col+kk,2) = limite((ii-1)*matrix_col,2)+(kk*L+x)/(v_design);
            else
                limite((ii-1)*matrix_col+kk,1) = limite((ii-1)*matrix_col+kk-1,2);
                limite((ii-1)*matrix_col+kk,2) = limite((ii-1)*matrix_col,2)+((kk+1)*L-x)/(v_design);
            end
            
        elseif ii==matrix_lin && ii==1% top/current layer which is 1st layer
            kk=kk+1;
            % odd filaments for even layers, even filaments for odd layers
            if mod(kk,2) == 1
                limite(kk,1) = (kk*L-x)/(v_design);
                limite(kk,2) = (kk*L+x)/(v_design);
                
            else
                limite(kk,1) = limite(kk-1,2);
                limite(kk,2) = ((kk+1)*L-x)/(v_design);
            end
        else%not 1st, and not last, i.e, layers in between
            kkk=kkk+1;
            % even filaments for even layers, odd filaments for odd layers
            if mod(kkk,2) == 1
                limite((ii-1)*matrix_col+kkk,1) = limite((ii-1)*matrix_col+kkk-1,2);
                limite((ii-1)*matrix_col+kkk,2) = limite((ii-1)*matrix_col,2)+((kkk)*L+x)/v(ii);
                
            % odd filaments for even layers, even filaments for odd layers
            else
                limite((ii-1)*matrix_col+kkk,1) = limite((ii-1)*matrix_col+kkk-1,2);
                limite((ii-1)*matrix_col+kkk,2) = limite((ii-1)*matrix_col,2)+((kkk+1)*L-x)/v(ii);
            end
        end
    end
end




% kk=0;kkk=0;
% % v is the optimum solution for the previous layer 
%                 % that is obtained a-prior
% for ii = 1:matrix_lin
%     for jj=1:matrix_col
%         if ii==matrix_lin % top/current layer 
%             kk=kk+1;
%             % odd filaments for even layers, even filaments for odd layers
%             if mod(kk,2) == 0
%                 if kkk==0
%                     limite(kkk+kk,1) = limite(kkk+kk-1,2);
%                     limite(kkk+kk,2) = (kk*L+x)/(v_design);
%                 else
%                     limite(kkk+kk,1) = limite(kkk+kk-1,2);
%                     limite(kkk+kk,2) = ((kkk+1)*L-x)/(v)+(2*kk)*(L-x)/(v_design);
%                 end
%                 
%                                 
%             % even filaments for even layers, odd filaments for odd layers
%             else % first enters here with numFilaments+1
%                 if kkk==0 && kk==1
%                     limite(kkk+kk,1) = (L-x)/(v_design);
%                     limite(kkk+kk,2) = (L+x)/(v_design);
%                 elseif kkk==0 && kk~=1
%                     limite(kkk+kk,1) = limite(kkk+kk-1,2);
%                     limite(kkk+kk,2) = (kk*L+x)/(v_design);
%                 else
%                     limite(kkk+kk,1) = limite(kkk+kk-1,2);
%                     limite(kkk+kk,2) = ((kkk+1)*L-x)/(v)+(2*kk)*(L-x)/(v_design);
%                 end
%                 
%             end
%         else
%             kkk=kkk+1;
%             % odd filaments
%             if mod(kkk,2) == 1
%                 limite(kkk,1) = (kkk*L-x)/(v);
%                 limite(kkk,2) = (kkk*L+x)/(v);
%             % even filaments
%             else
%                 limite(kkk,1) = limite(kkk-1,2);
%                 limite(kkk,2) = ((kkk+1)*L-x)/(v);
%             end
%         end
%     end
% end



% kk=0;kkk=0;
% for ii = 1:matrix_lin
%     for jj=1:matrix_col
%         if ii==matrix_lin % top/current layer 
%             kk=kk+1;
%             % odd filaments for even layers, even filaments for odd layers
%             if mod(kk,2) == 0
%                 limite(kkk+kk,1) = limite(kkk+kk-1,2);
%                 limite(kkk+kk,2) = ((kkk+1)*L-x)/(v)+(2*kk)*(L-x)/(v_design);
%             % even filaments for even layers, odd filaments for odd layers
%             else % first enters here with numFilaments+1
%                 if kkk==0 && kk==1
%                     limite(kkk+kk,1) = (L-x)/(v);
%                 else
%                     limite(kkk+kk,1) = limite(kkk+kk-1,2);
%                 end
%                 limite(kkk+kk,2) = ((kkk+1)*L-x)/(v)+(2*kk)*(L-x)/(v_design);
%             end
%         else
%             kkk=kkk+1;
%             % odd filaments
%             if mod(kkk,2) == 1
%                 limite(kkk,1) = (kkk*L-x)/(v);
%                 limite(kkk,2) = (kkk*L+x)/(v);
%             % even filaments
%             else
%                 limite(kkk,1) = limite(kkk-1,2);
%                 limite(kkk,2) = ((kkk+1)*L-x)/(v);
%             end
%         end
%     end
% end


% for kk = 1:1:(n+2)
%     if mod(kk,2) == 1
%         limite(kk,1) = (kk*L-x)/v;
%         limite(kk,2) = (kk*L+x)/v;
%     else
%         limite(kk,1) = limite(kk-1,2);
%         limite(kk,2) = ((kk+1)*L-x)/v;
%     end
% end

% temp=zeros(ceil((limite(n,2)+temp_extra)/dt)-1/dt,n); % BK
% for road = 1:n
%     line = 0;
%     for i = 0:dt:limite(n,2)
%         line = line + 1;
%         temp(line,road) = T_L;
%     end
% end



temp=zeros(ceil((limite(n,2)+temp_extra)/dt)-1/dt,n); % BK
for layerIdx = 1:matrix_lin
    for filamentIdx = 1:matrix_col
        line=0;
        if layerIdx==matrix_lin % top/current layer
            for i = 0:dt:limite(n,2)+temp_extra
                line=line+1;
                temp(line,(layerIdx-1)*matrix_col+filamentIdx) = T_design;
            end
        else
            for i = 0:dt:limite(n,2)+temp_extra
                line=line+1;
                % T_N is the optimum solution for the previous layer 
                % that is obtained a-prior
                temp(line,(layerIdx-1)*matrix_col+filamentIdx) = T_L(layerIdx);
            end
        end
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
                    if layer == matrix_lin
                        temp(p,num)=(T_design-Q(num))*exp(scalar*b(num)* ...
                            (t-limite(num,1))) + Q(num);
                    else
                        temp(p,num)=(T_L(layer)-Q(num))*exp(scalar*b(num)* ...
                            (t-limite(num,1))) + Q(num);
                    end
                end
                %Saving the last temperature of the period time of cooling down
                T_begin(num) = temp(p,num);
                %____________________________________ STEP 13 ____________________________________
%                 %Verification of the value of Biot Number
%                 Biot(num) = (vol/A_p)*(b(num)/kt(num));
%                 if Biot(num)>=0.1
%                     'WARNING! We cannot use a Lumped System';
%                 end
                
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
                            if layer == matrix_lin
                                save_T(num) = (T_design-Q(num))* ...
                                    exp(scalar*b(num)*(t-limite(num,1)))+ Q(num);
                            else
                                save_T(num) = (T_L(layer)-Q(num))* ...
                                    exp(scalar*b(num)*(t-limite(num,1)))+ Q(num);
                            end
                            
                            save_lim(1,1) = limite(num,1);
                            save_lim(1,2) = limite(num,1);
                        else
                            for j=1:num-2
                                save_T(j) = (T_begin(j)-Q(j))*exp(scalar*b(j)* ...
                                    (t-limite(j,1)))+Q(j);
                            end
                            save_T(num-1) = (T_begin(num-1)-Q(num-1))* ...
                                exp(scalar*b(num-1)*(t-limite(num,1)))+Q(num-1);
                            % BK
                            if layer == matrix_lin % if this is the top/last layer
                                save_T(num) = (T_design-Q(num))* ...
                                    exp(scalar*b(num)*(t-limite(num,1)))+ Q(num);
                            else
                                save_T(num) = (T_L(layer)-Q(num))* ...
                                    exp(scalar*b(num)*(t-limite(num,1)))+ Q(num);
                            end
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
%                 %____________________________________ STEP 18 ____________________________________
%                 %Verification of the Biot Number
%                 for j=1:num
%                     Biot(j) = (vol/A_p)*(b(j)/kt(j));
%                     if Biot(j)>=0.1
%                         'WARNING! We can not use a Lumped System';
%                         j
%                         Biot(j)
%                     end
%                 end
            end
        end
        % =================================================================== %
        % For the remaining layers -> BK
    else

        for num = sum(vector(1:layer-1,1))+1:sum(vector(1:layer,1))
            
            % =================================================================== %
            % if the filament is the first filament on that layer
            if num == sum(vector(1:layer-1,1))+1
                %____________________________________ STEP 9 ____________________________________
                a(num,2) = 1; %Activation of the contact down filament
                a(num-matrix_col,4) = 1; %Activation of the contact top filament
                
                % rest of the filaments on that layer
            else % if the filament is not the first filament on that layer
                
                %____________________________________ STEP 14 ____________________________________
                %Activation of the contacts (horizontal)
                a(num-1,3) = 1; a(num,1) = 1;
                %Activation of the contacts (vertical)
                a(num,2) = 1; a(num-matrix_col,4) = 1;
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
%             chk = size(save_T,2);
            
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
                            T(j,4) = save_T(j+matrix_col);
                            viz(4,j) = j+matrix_col;
                        end
                        if j > 1 && j < sum(vector(1,1))
                            T(j,1) = save_T(j-1);
                            viz(1,j) = j-1;
                            T(j,3) = save_T(j+1);
                            viz(3,j) = j+1;
                            T(j,4) = save_T(j+matrix_col);
                            viz(4,j) = j+matrix_col;
                        end
                        if j == sum(vector(1,1))
                            T(j,1) = save_T(j-1);
                            viz(1,j) = j-1;
                            T(j,4) = save_T(j+matrix_col);
                            viz(4,j) = j+matrix_col;
                        end

                        
                        % layers' first filament
                        if mod(j,matrix_col)==1 && j~=1
                            T(j,3) = save_T(j+1);
                            viz(3,j) = j+1;
                            T(j,2) = save_T(j-matrix_col);
                            viz(2,j) = j-matrix_col;
                            if layer-1 > 1
                                T(j-matrix_col,4) = save_T(j);
                                viz(4,j-matrix_col) = j;
                            end
                        end
                        % layers' inner filaments
                        if mod(j,matrix_col)~=0 && mod(j,matrix_col)~=1 && ...
                                j>matrix_col
                            T(j,1) = save_T(j-1);
                            viz(1,j) = j-1;
                            T(j,3) = save_T(j+1);
                            viz(3,j) = j+1;               
                            T(j,2) = save_T(j-matrix_col);
                            viz(2,j) = j-matrix_col;
                            if layer-1 > 1
                                T(j-matrix_col,4) = save_T(j);
                                viz(4,j-matrix_col) = j;
                            end
                        end
                        % layers' last filament
                        if mod(j,matrix_col)==0 && j~=matrix_col
                            T(j,1) = save_T(j-1);
                            viz(1,j) = j-1;        
                            T(j,2) = save_T(j-matrix_col);
                            viz(2,j) = j-matrix_col;
                            if layer-1 > 1
                                T(j-matrix_col,4) = save_T(j);
                                viz(4,j-matrix_col) = j;
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
                    % BK
                    if layer == matrix_lin % if this is the top/last layer
                        save_T(num) = (T_design-Q(num))*exp(scalar...
                            *b(num)*(t-limite(num,1)))+Q(num);
                    else
                        save_T(num) = (T_L(layer)-Q(num))*exp(scalar*...
                            b(num)*(t-limite(num,1)))+Q(num);
                    end
                    
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
            %                 %____________________________________ STEP 18 ____________________________________
            %                 %Verification of the Biot Number
            %                 for j=1:num
            %                     Biot(j) = (vol/A_p)*(b(j)/kt(j));
            %                     if Biot(j)>=0.1
            %                         'WARNING! We can not use a Lumped System';
            %                         j
            %                         Biot(j)
            %                     end
            %                 end
            
        end
    end
end


end
