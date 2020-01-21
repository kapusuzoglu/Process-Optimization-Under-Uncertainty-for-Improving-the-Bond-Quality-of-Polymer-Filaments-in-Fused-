% MCMC - Metropolis Hastings
%% actual measurement reading: the observed value


N_Samples = 1000;
% sigma_obs_initial = 5; % in mm ~ usually the bond radius is 0.05 mm

b2_ = 0.056; % Model parameter for temp. dependent viscosity
b2_l = 0.016; b2_r = 0.066;
roC = 1040*1290; % Model parameter for temp. dependent viscosity
roC_l = 1000*1100; roC_r = 1100*2000;
sigma_obs_initial = 0.05; % in mm ~ usually the bond radius is 0.05 mm
layer = 6;
% measured neck radii, written in interface order starting from the 1st
% layer in mm
% bond_test = [0.23368,0.37084,0.42672,0.49276,0.47244,0.51562,0.45466,...
%     0.52832,0.4064,0.4572,


Sample_Temp_Speed = 1.0*[234 34; 219 23; 242 26;
                         223 29; 247 37; 255 36;...
                         252 21];
                     
                     
% [t234s34_3_temps.mat; t219s23_4_temps.mat; t242s26_1_temps.mat; 
% t223s29_4_temps.mat; t247s37_1_temps.mat; t255s46_1_temps.mat;
% t252s21_1_temps.mat;]

bond_test = [0.38862,0.4318,0.42164,0.46736,0.47752,0.51816,0.49784,0.5207,0.51562,0.48514,0.4699,0.44704,0.40132,0.40132;...
0.44958,0.55372,0.56134,0.57404,0.5461,0.59182,0.53848,0.54356,0.49276,0.49276,0.46228,0.39624,0.39878,0.46228;...
0.3937,0.55626,0.56388,0.60706,0.56388,0.5334,0.4953,0.45466,0.4826,0.46736,0.4445,0.3175,0.35052,0.38608;...
0.39116,0.51562,0.52578,0.63246,0.60198,0.55118,0.47244,0.4826,0.46736,0.48006,0.46228,0.39116,0.38608,0.381;...
0.4318,0.5207,0.59436,0.50546,0.58928,0.53848,0.4826,0.4699,0.4572,0.4191,0.44958,0.42164,0.3683,0.39116;...
0.36576,0.42418,0.52578,0.52324,0.51816,0.53848,0.50292,0.49276,0.42418,0.4445,0.40132,0.39624,0.3556,0.34036;...
    

0.43173904,0.48947324,0.50493422,0.50115724,0.5785358,0.61870336,0.58698384,0.62997842,0.59888374,0.6147689,0.60409836,0.49989232,0.41950894,0.5190998;...
0.519049,0.5785739,0.58913014,0.58656982,0.56373014,0.56924956,0.56672988,0.60197746,0.60713874,0.61403738,0.53280564,0.56506872,0.45573442,0.42179748;...
0.46244256,0.61880242,0.58718196,0.6180201,0.61208158,0.61346588,0.58089038,0.58718196,0.5451729,0.56220868,0.5415026,0.5354193,0.4821936,0.4835652;...
0.462407,0.56802274,0.5750433,0.60108592,0.61307218,0.62429644,0.55603902,0.61515244,0.57609994,0.5856859,0.55785258,0.52633118,0.50043842,0.4701286;...
0.49094136,0.5489702,0.55971948,0.58943494,0.5979033,0.5921502,0.56146446,0.57645046,0.59953398,0.58283602,0.54763924,0.5282946,0.48722788,0.46245272;...
0.44813474,0.52808124,0.53608478,0.59736736,0.53349144,0.57643268,0.53484272,0.54762654,0.5466842,0.56845962,0.53246782,0.5312156,0.49306734,0.45228002;...

0.54150768,0.54807358,0.55617364,0.52196238,0.535178,0.56175656,0.6057011,0.63568072,0.63545974,0.66295778,0.5824982,0.58606944,0.51688238,0.54834028;...
0.56591708,0.6246876,0.62308486,0.64505586,0.61932566,0.62021974,0.61256164,0.60204096,0.62856872,0.69019166,0.67615816,0.65758568,0.58494422,0.59459622;...
0.54873398,0.62010544,0.6432677,0.67705732,0.69760592,0.66809874,0.65692782,0.67513708,0.6311773,0.64018668,0.64587628,0.6138799,0.59119008,0.60037472;...
0.57531,0.63385446,0.68160646,0.70896734,0.72745092,0.69035676,0.71760314,0.68299584,0.68746624,0.67827652,0.67199002,0.61706252,0.60833,0.60790836;...
0.52084986,0.67046602,0.6428105,0.65814194,0.62973712,0.64737996,0.68527422,0.71647304,0.6893941,0.61033406,0.65337944,0.69042534,0.64621156,0.62503304;...
0.49080928,0.59097164,0.58768996,0.61816234,0.64793622,0.64242188,0.62777624,0.66053462,0.6208141,0.64635888,0.68806822,0.6163056,0.58142886,0.58870596;...

0.4350004,0.500507,0.5303901,0.5480177,0.56173878,0.5689981,0.56482742,0.55253636,0.5595874,0.58387234,0.5161915,0.4612513,0.43331384,0.42196766;...
0.44862242,0.58602118,0.61177678,0.61554868,0.61158628,0.60180982,0.57809638,0.57933844,0.52659534,0.54688994,0.53090064,0.47236634,0.44041314,0.4378452;...
0.41031922,0.51967384,0.54105556,0.58049414,0.58480452,0.61651642,0.55451756,0.50642266,0.4874387,0.527685,0.45568362,0.42885868,0.43361102,0.44055284;...
0.3909949,0.52229004,0.56043068,0.53613558,0.56969914,0.56820816,0.51941222,0.51244246,0.47563786,0.48608742,0.42275252,0.3707257,0.35813238,0.37316664;...
0.4121785,0.5350637,0.590423,0.56943752,0.59579002,0.52862988,0.4997196,0.47566834,0.50154332,0.5381117,0.47614586,0.35383216,0.30835092,0.34651696;...
0.38857936,0.49550828,0.54940708,0.56415686,0.53723032,0.54927246,0.51648614,0.49673256,0.48010064,0.50453798,0.47772828,0.44153836,0.3950335,0.34187638;...


0.48697642,0.54173628,0.58204354,0.6334506,0.71483728,0.70277228,0.686181,0.685673,0.69745606,0.7117207,0.6978523,0.65822068,0.6735318,0.64583056;...
0.44146216,0.5521833,0.56614314,0.63512192,0.6039866,0.55593488,0.55131462,0.57001918,0.55765192,0.56694324,0.5634482,0.5281803,0.51191414,0.520954;...
0.44977304,0.49388776,0.5610225,0.60408566,0.56909462,0.59820302,0.57758838,0.54587394,0.56660796,0.53937408,0.49987708,0.47286418,0.43899582,0.46756828;...
0.47943516,0.54695598,0.60337192,0.61815472,0.63889128,0.59216798,0.5971667,0.54893972,0.56864758,0.60009786,0.51827938,0.47008034,0.45374814,0.46762416;...
0.44791884,0.56211724, 0.5443347,0.57093612,0.4236974,0.59906408,0.48954944,0.50017172,0.47689516,0.54124352,0.45867574,0.47834296,0.37190172,0.44401232;...
0.4414266,0.51579272,0.54539642,0.6175502,0.59111134,0.5808218,0.52012596,0.51555904,0.55867808,0.55106316,0.49651158,0.502412,0.4172331,0.37376862;...

0.4003802,0.46118272,0.52152296,0.55893716,0.5373243,0.54146196,0.52077874,0.51377596,0.50603658,0.4873244,0.4704969,0.40757856,0.33943544,0.34058098;...
0.45451522,0.52318666,0.55244492,0.59022996,0.54899306,0.50519076,0.31620206,0.38583108,0.38320472,0.42009822,0.38410388,0.34455354,0.28713684,0.28238196;...
0.45597318,0.59106816,0.5774944,0.62265306,0.57917334,0.47201328,0.25991312,0.3564382,0.3927094,0.37080952,0.36638992,0.34158174,0.30037786,0.25989026;...
0.42583608,0.5380609,0.58333386,0.51698144,0.6286373,0.46335188,0.35287204,0.36197794,0.32757618,0.39344854,0.36700714,0.37292026,0.26477722,0.23824692;...
0.4367911,0.45451522,0.52296314,0.59928252,0.60136024,0.44869608,0.39428166,0.33972246,0.35554666,0.4014597,0.41152826,0.34999168,0.26832306,0.24936704;...
0.39101268,0.44653708,0.51600608,0.50715418,0.57846976,0.51341274,0.34140902,0.40498014,0.38158166,0.44459652,0.4210939,0.37876988,0.27357578,0.26587704;...


0.5518912,0.6269736,0.657606,0.68438776,0.6373876,0.6759194,0.70395338,0.61932312,0.70766178,0.70857872,0.64976502,0.6629146,0.62848236,0.74328782;...
0.57904888,0.59551824,0.6766306,0.66521838,0.68699888,0.72004174,0.78422754,0.72168004,0.69519546,0.69670168,0.59077352,0.63762128,0.53293772,0.71866506;...
0.56349646,0.63391034,0.66410586,0.65397634,0.64886078,0.67883532,0.67612006,0.64920368,0.62444884,0.64120776,0.63840106,0.63662306,0.63966852,0.66999104;...
0.5063109,0.61478668,0.71725028,0.69045836,0.76191364,0.64057784,0.73488042,0.70996302,0.70214744,0.65829688,0.6385814,0.62281054,0.581787,0.60792614;...
0.55894986,0.66879978,0.44379642,0.70409054,0.72319896,0.75304904,0.70470268,0.71902828,0.66672968,0.65223898,0.58577226,0.56827166,0.59377834,0.61651896;...
0.48673004,0.58683906,0.66327528,0.63945262,0.6787769,0.66257424,0.72108822,0.69915278,0.75986894,0.73452736,0.69030596,0.61689742,0.58194448,0.65120012]/2;



% numInt = size(bond_test,2);

% 
% % import temperature data
% Temp.S1 = load('.\Experiment_Analytical_Abaqus\T234S34-selected\t234s34_3_temps.mat');
% Temp.S2 = load('.\Experiment_Analytical_Abaqus\T219S23-selected\t219s23_4_temps.mat');
% Temp.S3 = load('.\Experiment_Analytical_Abaqus\T242S26-selected\t242s26_1_temps.mat');
% Temp.S4 = load('.\Experiment_Analytical_Abaqus\T223S29-selected\t223s29_4_temps.mat');
% Temp.S5 = load('.\Experiment_Analytical_Abaqus\T247S37-selected\t247s37_1_temps.mat');
% Temp.S6 = load('.\Experiment_Analytical_Abaqus\T255S46-selected\t255s46_1_temps.mat');
% Temp.S7 = load('.\Experiment_Analytical_Abaqus\T252S21-selected\t252s21_1_temps.mat');



NumTest = size(bond_test,1)/layer;

%%
tic
[MC_samples, AcccRate] = Metro(N_Samples, b2_, ...
                    bond_test,numInt, b2_l, b2_r, NumTest,...
                    roC,roC_l,roC_r,Sample_Temp_Speed);
toc
%%



% ======================================================================= %
function [MC_chains, accept_rate] = Metro(Nsamples, b2_initial, ....
                bond_test, n_test, b2_left, b2_right, NumTest,...
                    roC,roC_l,roC_r,Sample_Temp_Speed)
% Metropolis for posterior samples of C and sigma_eps

% INPUTS
% Nsamples: total number of samples to draw
% b2_inital: parameter in Neck growth Model, initial
% sigma_eps: eps = bond_obs - bond_pred, eps ~ N(0, sigma_eps), passed in
% bond_test, n_test: bond radius measurement from the lab test
% b2_left, b2_right: parameters for the prior of b2

% OUTPUTS
% MC_chains(s,:): the sâ€™th sample (of size d)
% accept_rate: ratio of accepted moves

% place holder
MC_chains = zeros(Nsamples, 2);
x = [b2_initial,roC];
naccept = 0;

% evaluate the target pdf value at the initial point
logpOld = LN_TargetPDF(b2_initial, bond_test,n_test, ...
                        b2_left, b2_right, NumTest,...
                    roC,roC_l,roC_r,Sample_Temp_Speed);

for t = 1: Nsamples
    
    % propose a new point - Metropolis Hastings algorithm
    
    % to ensure C, sigma_eps are positive values, use a lognormal proposal,
    % assume the proposal lognormal has a mean of current point, sd 0.5*
    % calculate the parameters for lognormal distribution
    
    b2_curr = x(1, 1);
    % mu_C = log(C_curr^2 / sqrt((0.5 * C_curr)^2 + C_curr^2));
    mu_b2 = log(b2_curr) - 0.5 * log(1.25);
    % sig_C = sqrt(log((0.5 * C_curr)^2 / C_curr^2 + 1));
    sig_b2 = sqrt(log(1.25));
    
    b2_p = lognrnd(mu_b2, sig_b2);
    
    
    roc_curr = x(1, 2);
    % mu_C = log(C_curr^2 / sqrt((0.5 * C_curr)^2 + C_curr^2));
    mu_roc = log(roc_curr) - 0.5 * log(1.25);
    % sig_C = sqrt(log((0.5 * C_curr)^2 / C_curr^2 + 1));
    sig_roc = sqrt(log(1.25));
    
    roc_p = lognrnd(mu_roc, sig_roc);
    
%     sigma_eps_curr = x(1, 2);
%     mu_sigma = log(sigma_eps_curr) - 0.5 * log(1.25);
%     sig_sigma = sqrt(log(1.25));
%     
%     sigma_eps_p = lognrnd(mu_sigma, sig_sigma);
    
    % evaluate the target pdf value at the proposed point
    logpProposed = LN_TargetPDF(b2_p, bond_test,n_test, ...
                        b2_left, b2_right, NumTest,...
                    roc_p,roC_l,roC_r,Sample_Temp_Speed);
    
    % Proposal at current point                
    p1 = lognpdf(b2_curr, mu_b2, sig_b2)* ...
         lognpdf(roc_p, mu_roc, sig_roc);
    % Proposal at proposed point
    p2 = lognpdf(b2_p, mu_b2, sig_b2) * ...
         lognpdf(roc_p, mu_roc, sig_roc);
    % for mu_epsF and sigma_epsF, since the proposal pdf is uniform, it's
    % the same at both current point and the proposed point, no need to
    % include in p1 p2.
    
    % acceptance ratio
    alpha = exp(logpProposed - logpOld) * p1/p2;
    
    r = min(1, alpha);
    u = rand(1,1);
    if u < r % accept the proposed point as the new point
        x = [b2_p, roc_p];
        naccept = naccept + 1;
        logpOld = logpProposed;
    end
    MC_chains(t, :) = x;
    if mod(t, 500) == 0
        filename1 = ['.\MCMC_TempModelandBondModel\MC_', num2str(t), 'Samples.txt'];
        dlmwrite(filename1, MC_chains(1:t, :));
        filename2 = ['.\MCMC_TempModelandBondModel\MC_', num2str(t), 'Samples_NumofAcceptedPts.txt'];
        dlmwrite(filename2, naccept);
    end
end
accept_rate = naccept / Nsamples;
end
%%

% ======================================================================= %
% Natural log value of the target pdf evaluated at b2 and sigma_eps

% INPUTS
% b2_inital: parameter in Neck growth Model, initial
% sigma_eps: eps = bond_obs - bond_pred, eps ~ N(0, sigma_eps), passed in
% bond_test, n_test: bond radius measurement from the lab test
% b2_left, b2_right: parameters for the prior of b2

% OUTPUT
% LN_PDF: the natural log value of the target pdf evaluated at b2 and sigma_eps
%    the target pdf is propotional to the posterior
function LN_PDF = LN_TargetPDF( b2_, bond_test,n_int, ...
                            b2_left, b2_right, NumTest,...
                    roC,roC_l,roC_r,Sample_Temp_Speed)

% ======================================================================= %
%__________________________       INPUTS     ______________________________
%Process Variables
%__________________________________________________________________________
%Filament dimensions
%Process Variables
%________________________________________________________________________
% T_L = 215; %Extrusion temperature (ºC)
T_E = 110; %Temperature of the envelope (ºC)
% v = 0.032; %Velocity of the extrusion head (m/sec)
w = 0.8e-3; % layer width [meters]
h = 0.7e-3; % layer height [meters]
L = 0.035; %Length of the filament (meters)

numLayers = 6;
num_Filaments = n_int+1;
% numInt = numFilaments-1;


% Process parameters
% ProcessParam = [wth, hth, layer, num_Filaments];
ProcessParam = [T_E, w, h, L, num_Filaments,numLayers];

%__________________________________________________________________________
% Material parameters for the bond length model
Surf_tens = 0.29; % Surface tension
b1_ = 0.00345; % Model parameter for temp. dependent surface tension
Visco = 5100; % viscosity
% b2_ = 0.056; % Model parameter for temp. dependent viscosity
BondModelParam = [Surf_tens, b1_, Visco, b2_];
%__________________________________________________________________________


% The low fidelity 1-D analytical model
% bond_radii: bond length between two filaments:
%        Output of the heat transfer model that does not consider axial
%        heat conduction & neck growth model, in mm.
[bond_radii] = Neck_Growth_Model(ProcessParam, roC, BondModelParam, NumTest,Sample_Temp_Speed);


% ------ epsilon: the difference between prediction and measurements ------
epsilon = bond_test - bond_radii;

% ------ Likelihood ------
likelihood = normpdf(epsilon, 0, 0.05);

% ------ Priors ------
% Prior of b2: uniform, pdf_b2 = 1/(b2_right - b2_left).
% Prior of sigma_eps: (uninformed, Jeffreys Prior, propotional to 1/sigma),
%                     density = 1/sigma_eps.
%                     This is an IMPROPER prior.

% ------ logPropTarget = log(Likelihoods * Priors)
%                      = sum[log(Likelihood)] + sum[log(prior)] ------
% All measurements are independent, 
% LN_PDF = sum(log(likelihood)) + log(1/(b2_right - b2_left)) + log(1/sigma_eps)
LN_PDF = sum(log(likelihood)) - log(b2_right - b2_left) - log(roC_r - roC_l);

end





% ======================================================================= %

function[AllBondLengths] = Neck_Growth_Model(ProcessParam,RoC, inputs, NumTest,Sample_Temp_Speed)

%Filament dimensions
w = ProcessParam(2); % Layer Thickness (meters)
h = ProcessParam(3); % Layer Height (meters)
L = ProcessParam(4); % Layer Length (meters)
numFilaments = ProcessParam(5); % considered interface number on that layer
layer = ProcessParam(6); % considered layer number
% numFilaments - 1 = interface number (ranging from 1-14 for 15 filaments)
% for that given layer

% convert Celcius to Kelvin
KK = 273.15; % 0 Celcius in K



% ======================================================================= %
%% Temperature Model
%% This will plot the temperatures at line B (x = 175 pixels)
Matrix = ones(layer,numFilaments);
x = L/2; % location of the cut in meters



HalfBondLengths=zeros(layer,numFilaments-1);
AllBondLengths=zeros(NumTest*layer,numFilaments-1);
NeckRadius=zeros(1,numFilaments-1);

for iii=1:NumTest
    % nozzle temperature and printer speed for each sample
    TV=Sample_Temp_Speed(iii,:);
    % the nozzle temperature is observed to be at least 10 degrees less
    T_N = TV(1)-10; 
    v_p = TV(2);
    [Temp_C, time] = FDM_TempModel(Matrix,x,T_N,v_p,ProcessParam,RoC);
    
    % convert Celcius to Kelvin
    Temp = Temp_C + KK;
    
    % Average two adjacent lines' temperatures to obtain interface temperatures
    % for each interface on that layer
    int_Temp = zeros(size(Temp,1),numFilaments-1);
    T_idx = zeros(1,numFilaments-1);
    T_idx_ = zeros(1,numFilaments-1);
    for ii=2:numFilaments
        T_idx(ii-1) = find(Temp(:,ii)<T_N+KK, 1, 'first');
        % find the first index that temp value is lower thatn 145 celcius
        T_idx_(ii-1) = find(Temp(:,ii)<145+KK, 1, 'first');
        int_Temp(T_idx(ii-1):T_idx_(ii-1),ii-1) = (Temp(T_idx(ii-1):T_idx_(ii-1),ii-1)+Temp(T_idx(ii-1):T_idx_(ii-1),ii))/2;
    end
    
    int_Temp_model = int_Temp'; % interface temperature

    tfinal=time(end);
    num=size(int_Temp_model,1)-T_idx(1)+1;
    dt=tfinal/num;

        
    for jj=1:layer
        
        for rr=1:numFilaments-1 % loop over each interface for layer 1

            % interface temperature for experiment
            T = int_Temp_model(rr,:);
            
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
            
            % ao = w/2; % initial radius: half of layer width in meters
            %     ao = sqrt(w/2*h/2); % initial radius: half of layer width in meters
            aa=w/2;bb=h/2;
            ao = aa*bb/sqrt(aa^2*sin(0.25*pi)^2+bb^2*cos(0.25*pi)^2);
            
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
            bondLength=ao*sin(theta);
            
            % find the time when temperature reaches 112 celcius
            index1=find(Ref_Temp+KK>T,1,'first');
            [~,index2]=min(abs(time(index1)-t_bond));
            
%             if rr==1
%                 % find the time index that corresponds to the same time
%                 [~,index2]=min(abs(time(index1)-t_bond));
%             else
%                 % find the time index that corresponds to the same time
%                 [~,index2]=min(abs(exp_time(index1+firstEl-1)-t_bond));
%             end
            
            % save neck radius for each interface on a given layer
            NeckRadius(rr) = bondLength(index2)*1e3; % neck radius in mm
        end
    
        HalfBondLengths(jj,:) = NeckRadius;
    end
    AllBondLengths((iii-1)*layer+1:iii*layer,:) = HalfBondLengths;
end


end


    
function[temp,abcissa] = FDM_TempModel(matrix,x,T_n,V_p,ProcessParam,MaterialParam)

%Process Variables
T_L = T_n; %Extrusion temperature (ºC)
T_E = ProcessParam(1); %Temperature of the envelope (ºC)
v = V_p; % printer speed/Velocity of the extrusion head (m/sec)
%Filament dimensions
Wid = ProcessParam(2); % Layer Thickness (meters)
Ht = ProcessParam(3); % Layer Height (meters)
L = ProcessParam(4); % Length of the filament (meters)

roC = MaterialParam; % Density (kg/m^3)
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
h(1,1) = 86; % filament and left adjacent filament
h(1,2) = 86; % filament and down adjacent filament
h(1,3) = 86; % filament and right adjacent filament
h(1,4) = 86; % filament and top adjacent filament
h(1,5) = 86; % filament and support
% h(1,5) = h_cond;
%Fraction of perimeter contact between
lambda(1,1) = 0.1; % filament and left adjacent filament
lambda(1,2) = 0.1; % filament and down adjacent filament
lambda(1,3) = 0.2; % filament and right adjacent filament
lambda(1,4) = 0.1; % filament and top adjacent filament
lambda(1,5) = 0.1; % filament and support
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
    
    