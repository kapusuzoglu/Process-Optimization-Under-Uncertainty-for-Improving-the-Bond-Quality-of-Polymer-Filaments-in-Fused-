% MCMC - Metropolis Hastings
%% actual measurement reading: the observed value

% read the csv file
M = csvread('032_Dec_7.csv');
% diameter of the intial hole + introduced crack, deemed as a0 at n0 = 0.
a0 = M(1, 2);

% n means Cycles 
n_test = M(:, 1);

% a_raw means Crack Length (excludes the diameter of the intial hole)
a_raw = M(2:end, 2);
a_test = zeros(size(n_test));
a_test(1, 1) = a0;
a_test(2:end, 1) = a_raw + a0;

% start at a long crack
a_start = a_test(2, 1);
% corresponding cycle count
n_start = n_test(2, 1);
% end poit of test
nf = n_test(end);

%% Load the GP model for K1
load('gprMdl_K1.mat');

%%
tic
[MC_samples, AcccRate] = ...
    Metro(10, ...
          3.4e-9, 0.01, 3e-9, 9e-10, gprMdl_K1, 67, 3.5241,...
          5000, 2500, a_start, n_start, nf, ...
          a_test, n_test,...
          1.8e-9, 5.4e-9, 2e-9, 4e-9);
toc
%%



% ======================================================================= %
function [MC_chains, accept_rate] = ...
         Metro(Nsamples, ...
               C_initial, sigma_eps_initial, mu_epsF_initital, ...
               sigma_epsF_initial, gprMdl_K1, Kc, m, ...
               Fmax_lbs, Fmin_lbs, a0, n0, nf, ...
               a_test, n_test, ...
               C_left, C_right, mu_epsF_right, mu_epsF_left)
% Metropolis for posterior samples of C and sigma_eps

% INPUTS
% Nsamples: total number of samples to draw
% C_inital: parameter in Forman's Model, initial
% sigma_eps_initial: eps = a_obs - a_pred, eps ~ N(0, sigma_eps), initial
% mu_epsF_initital, sigma_epsF_initial:
%    epsF ~ N(mu_epsF, sigma_epsF), these are inital values of that.
% gprMdl_K1: GP model for stress intensity factor K1
% Kc: parameter in Forman's model, for AL7075_T6, use 67, MPa*sqrt(m)
% m: parameter in Forman's model, for AL7075_T6, use 3.5241
% Fmax_lbs, Fmin_lbs: max min loads of the constant amplitude loding, lbs
% a0: starting point of the prediction of crack growth path
% n0: starting cycle count
% nf: ending cycle count of prediction
% a_test, n_test: crack size measurement a at cycle n from the lab test
% C_left, C_right: parameters for the prior of C

% OUTPUTS
% MC_chains(s,:): the sâ€™th sample (of size d)
% accept_rate: ratio of accepted moves

% place holder
MC_chains = zeros(Nsamples, 4);
x = [C_initial, sigma_eps_initial, mu_epsF_initital, sigma_epsF_initial];
naccept = 0;

% evaluate the target pdf value at the initial point
logpOld = logp(C_initial, sigma_eps_initial, mu_epsF_initital,...
              sigma_epsF_initial, gprMdl_K1, Kc, m, ...
              Fmax_lbs, Fmin_lbs, a0, n0, nf, ...
              a_test, n_test, ...
              C_left, C_right, mu_epsF_right, mu_epsF_left);

for t = 1: Nsamples
    
    % propose a new point - Metropolis Hastings algorithm
    
    % to ensure C, sigma_eps are positive values, use a lognormal proposal,
    % assume the proposal lognormal has a mean of current point, sd 0.5*
    % calculate the parameters for lognormal distribution
    
    C_curr = x(1, 1);
    % mu_C = log(C_curr^2 / sqrt((0.5 * C_curr)^2 + C_curr^2));
    mu_C = log(C_curr) - 0.5 * log(1.25);
    % sig_C = sqrt(log((0.5 * C_curr)^2 / C_curr^2 + 1));
    sig_C = sqrt(log(1.25));
    
    C_p = lognrnd(mu_C, sig_C);
    
    sigma_eps_curr = x(1, 2);
    mu_sigma = log(sigma_eps_curr) - 0.5 * log(1.25);
    sig_sigma = sqrt(log(1.25));
    
    sigma_eps_p = lognrnd(mu_sigma, sig_sigma);
    
    % proposal of mu and sigma of epsF are to be drawn from a uniform
    % distribution: centered at the current point, left end at 0.
    mu_epsF_curr = x(1, 3);
    mu_epsF_p = rand * mu_epsF_curr * 2;
    
    sigma_epsF_curr = x(1, 4);
    sigma_epsF_p = rand * sigma_epsF_curr * 2;
        
    % evaluate the target pdf value at the proposed point
    logpProposed = LN_TargetPDF(C_p, sigma_eps_p, mu_epsF_p, sigma_epsF_p,...
                        gprMdl_K1, Kc, m, ...
                        Fmax_lbs, Fmin_lbs, a0, n0, nf, ...
                        a_test, n_test, ...
                        C_left, C_right, mu_epsF_right, mu_epsF_left);
    
    % Proposal at current point                
    p1 = lognpdf(C_curr, mu_C, sig_C) * ...
         lognpdf(sigma_eps_curr, mu_sigma, sig_sigma);   
    % Proposal at proposed point
    p2 = lognpdf(C_p, mu_C, sig_C) * ...
         lognpdf(sigma_eps_p, mu_sigma, sig_sigma);
    % for mu_epsF and sigma_epsF, since the proposal pdf is uniform, it's
    % the same at both current point and the proposed point, no need to
    % include in p1 p2.
    
    % acceptance ratio
    alpha = exp(logpProposed - logpOld) * p1/p2;
    
    r = min(1, alpha);
    u = rand(1,1);
    if u < r % accept the proposed point as the new point
        x = [C_p, sigma_eps_p, mu_epsF_p, sigma_epsF_p];
        naccept = naccept + 1;
        logpOld = logpProposed;
    end
    MC_chains(t, :) = x;
    if mod(t, 5) == 0
        filename1 = ['MC_', num2str(t), 'Samples.txt'];
        dlmwrite(filename1, MC_chains(1:t, :));
        filename2 = ['MC_', num2str(t), 'Samples_NumofAcceptedPts.txt'];
        dlmwrite(filename2, naccept);
    end
end
accept_rate = naccept / Nsamples;
end
%%

% ======================================================================= %
% Natural log value of the target pdf evaluated at C and sigma_eps

% INPUTS
% C: parameter in Forman's Model, passed in
% sigma_eps: eps = a_obs - a_pred, eps ~ N(0, sigma_eps), passed in
% gprMdl_K1: GP model for stress intensity factor K1
% Kc: parameter in Forman's model, for AL7075_T6, use 67, MPa*sqrt(m)
% m: parameter in Forman's model, for AL7075_T6, use 3.5241
% Fmax_lbs, Fmin_lbs: max min loads of the constant amplitude loding, lbs
% a0: starting point of the prediction of crack growth path
% n0: starting cycle count
% nf: ending cycle count of prediction
% a_test, n_test: crack size measurement a at cycle n from the lab test
% C_left, C_right: parameters for the prior of C

% OUTPUT
% LN_PDF: the natural log value of the target pdf evaluated at C and sigma_eps
%    the target pdf is propotional to the posterior
function LN_PDF = LN_TargetPDF(b2_, sigma_eps, b2_left, b2_right
                            bond_test, n_test)
              

% ======================================================================= %
%% Costa Model - Analytical
%__________________________       INPUTS     ______________________________
%Process Variables
%__________________________________________________________________________
T_ext = 215; %Extrusion temperature (ºC)
T_env = 110; %Temperature of the envelope (ºC)
v_printer = 0.032; %Velocity of the extrusion head (m/sec)
%Filament dimensions
wth = 0.8e-3; % layer width [meters]
hth = 0.7e-3; % layer height [meters]
Len = 0.035; %Length of the filament (meters)

numLayers = 6;
number_filament = 15; 
Matrix = ones(numLayers,number_filament);
x = Len/2; % location of the cut in meters

% Process parameters
ProcessParam = [T_ext, T_env, v_printer, wth, hth, Len];
%__________________________________________________________________________
% Material properties - ABS polymer
ro_ = 1040; % Density (kg/m^3)
C_ = 1290; % Specific heat (J/kg.K)
roC_ = ro_*C_; % Density*Specific heat
h_conv_ = 86; % Heat transfer coefficient (loss of heat by natural convection)

%Fraction of perimeter contact between
lam1_ = .2; % filament and left adjacent filament
lam2_ = .2; % filament and down adjacent filament
lam3_ = .2; % filament and right adjacent filament
lam4_ = .2; % filament and top adjacent filament
lam5_ = .2; % filament and support
lam_ = [lam1_ lam2_ lam3_ lam4_ lam5_];

h_conductivity = 0.15; % conductivity coefficient

%Thermal contact conductances between
h1 = 86; % filament and left adjacent filament
h2 = 86; % filament and down adjacent filament
h3 = 86; % filament and right adjacent filament
h4 = 86; % filament and top adjacent filament
h5 = 86; % filament and support
h_tc = [h1 h2 h3 h4 h5];
%__________________________________________________________________________
% Material parameters for the bond length model
Surf_tens = 0.29; % Surface tension
b1_ = 0.00345; % Model parameter for temp. dependent surface tension
Visco = 5100; % viscosity
% b2_ = 0.056; % Model parameter for temp. dependent viscosity
BondModelParam = [Surf_tens, b1_, Visco, b2_];
%__________________________________________________________________________
% Material parameters 
MaterialParam = [roC_, h_conv_, lam_, h_conductivity, h_tc, BondModelParam];
%__________________________________________________________________________


% The low fidelity 1-D analytical model
% bond1: bond length between two filaments:
%        Output of the heat transfer model that does not consider axial
%        heat conduction & neck growth model, in mm.
[bond1] = FDM_TempModel(Matrix, x, ProcessParam, MaterialParam);


% ------ epsilon: the difference between prediction and measurements ------
epsilon = bond_test(2:end) - bond1(n_test(2:end) - n0 + 1);

% ------ Likelihood ------
likelihood = normpdf(epsilon, 0, sigma_eps);

% ------ Priors ------
% Prior of C: uniform, pdf_C = 1/(C_right - C_left).
% Prior of sigma_eps: (uninformed, Jeffreys Prior, propotional to 1/sigma),
%                     density = 1/sigma_eps.
%                     This is an IMPROPER prior.
% Prior of mu_epsF: uniform, pdf = 1/(mu_epsF_right - mu_epsF_left)
% Prior of sigma_epsF: Jeffreys Prior, improper, density =  1/sigma_epsF

% ------ logPropTarget = log(Likelihoods * Priors)
%                      = sum[log(Likelihood)] + sum[log(prior)] ------
% All measurements are independent, 
% LN_PDF = sum(log(likelihood)) + log(1/(C_right - C_left)) + log(1/sigma_eps) +
%      log(1/(mu_epsF_right - mu_epsF_left)) + log(1/sigma_epsF)
LN_PDF = sum(log(likelihood)) - log(b2_right - b2_left) - log(sigma_eps);

end





% ======================================================================= %

function[bondLength] = FDM_TempModel(matrix,x,ProcessParam,MaterialParam)

%Process Variables
T_L = ProcessParam(1); %Extrusion temperature (ºC)
T_E = ProcessParam(2); %Temperature of the envelope (ºC)
v = ProcessParam(3); % printer speed/Velocity of the extrusion head (m/sec)
%Filament dimensions
Wid = ProcessParam(4); % Layer Thickness (meters)
Ht = ProcessParam(5); % Layer Height (meters)
L = ProcessParam(6); % Length of the filament (meters)

                
%____________________________________ STEP 1 ____________________________________
%Definition of the vector that contains the number of total filaments in each layer
matrix_lin = size(matrix,1);
matrix_col = size(matrix,2);
vector = zeros(matrix_lin,2);
ctr = 0; % # of layers
for ii = matrix_lin:-1:1
    ctr = ctr + 1;
    for j = 1:matrix_col
        if matrix(ii,j) ~= 0
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
temp_extra = 1; %Additional time computation after construction of the part
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


%____________________________________ STEP 5 ____________________________________
% Material Properties
roC = MaterialParam(1); % Density*Specific Heat (kg/m^3)*(J/kg.K)
h_conv = MaterialParam(2); % Heat transfer coefficient (loss of heat by natural convection)
lambda(1,1) = MaterialParam(3); % filament and left adjacent filament
lambda(1,2) = MaterialParam(4); % filament and down adjacent filament
lambda(1,3) = MaterialParam(5); % filament and right adjacent filament
lambda(1,4) = MaterialParam(6); % filament and top adjacent filament
lambda(1,5) = MaterialParam(7); % filament and support
conductivity = MaterialParam(8); % %Thermal conductivity (W/m.K)

%Thermal contact conductances between
h(1,1) = MaterialParam(9); % filament and left adjacent filament
h(1,2) = MaterialParam(10); % filament and down adjacent filament
h(1,3) = MaterialParam(11); % filament and right adjacent filament
h(1,4) = MaterialParam(12); % filament and top adjacent filament
h(1,5) = MaterialParam(13); % filament and support


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
% ctr = 0;
% number_filament = 0;
% for i = matrix_lin:-1:1
%     ctr = ctr + 1;
%     % mod(N,2)==0 % iseven
%     % mod(N,2)==1 % isodd
%     %     if isodd(contar) == 1
%     if mod(ctr,2) == 1
%         for j = 1:matrix_col
%             if matrix(i,j) ~= 0
%                 number_filament = number_filament + 1;
%                 scalar(number_filament) = -per/(ro(matrix(i,j))*area*C(matrix(i,j)));
%                 esc(number_filament) = h_conv/(ro(matrix(i,j))*L*C(matrix(i,j)));
%                 kt(number_filament) = conductivity(matrix(i,j));
%             end
%         end
%     else
%         for j = matrix_col:-1:1
%             if matrix(i,j) ~= 0
%                 number_filament = number_filament + 1;
%                 scalar(number_filament) = -per/(ro(matrix(i,j))*area*C(matrix(i,j)));
%                 esc(number_filament) = h_conv/(ro(matrix(i,j))*L*C(matrix(i,j)));
%                 kt(number_filament) = conductivity(matrix(i,j));
%             end
%         end
%     end
% end


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
    for ii = 0:dt:limite(n,2)
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



KK = 273.15; % 0 Celcius in K
% convert Celcius to Kelvin
Temp = temp + KK;

% 1st interface: numFilaments-1 gives the interface of concern
numFilaments = 2; 


int_Temp = zeros(size(Temp,1),numFilaments-1);
T_idx = zeros(1,numFilaments-1);
for ii=2:numFilaments
    % find the first index that represents the start of the extrusion
    T_idx(ii-1) = find(Temp(:,ii)<T_L+KK, 1, 'first');
    LastInd = find(Temp(T_idx(ii-1):end,ii)==T_L+KK, 1, 'first');
    % Average two adjacent lines' temperatures to obtain interface temperatures
    % for the 1st interface on that layer
    int_Temp(T_idx(ii-1):LastInd-1,ii-1) = (Temp(T_idx(ii-1):LastInd-1,ii-1)+Temp(T_idx(ii-1):LastInd-1,ii))/2;
end

t_final=abcissa(end);

% compute the time step
rr=numFilaments-1; % # of interface
num=size(int_Temp,1)-T_idx(rr)+1;
dt=t_final/num;

% Calculate the bond length at that corresponding interface
[bondLength]=Neck_Growth_Model(int_Temp(T_idx(rr):end,rr),dt,t_final,Wid,Ht,MaterialParam(14:end));



    function[NeckRadius]=Neck_Growth_Model(T,dt,tfinal,w,h,inputs)

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

    % ao = w/2; % initial radius: half of layer width in meters
    ao = sqrt(w/2*h/2); % initial radius: half of layer width in meters

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
    
    Kel = 273.15; % Kelvin conversion
    [~,idx]=min(abs(T-(240+Kel)));
    if T(1)>240+Kel
        Delta_Gamma(1:idx) = -0.0005;
    % Delta_Gamma(1:idx) = -0.000000001;
    end

    T_r_C = 240; % reference temperature in C
    T_r = T_r_C+Kel; % in K

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
    NeckRadius=ao*sin(theta);

    end




end


