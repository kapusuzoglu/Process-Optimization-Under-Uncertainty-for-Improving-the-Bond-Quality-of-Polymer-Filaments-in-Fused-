%==========================================================================    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%   Copyright: Berkcan Kapusuzoglu %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Heat transfer model + Bond length model.
%==========================================================================    
clc; clear; close all;

rng default % For reproducibility




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

%__________________________________________________________________________
% For Heat Transfer Model
L = 0.035; % filament lengths in meters
x_cut = L/2; % location of the cut in meters
T_bp = 110; % build plate temperature
% Process parameters
ProcessParam = [wth, hth, layer, numFilaments, T_bp, L];
%__________________________________________________________________________

format long

% % Load experimental temperature data for that sample
% temps = load('..\Experiment_Analytical_Abaqus\T240S42\t240s42_1_temps.mat');

% % Measured Bond Lengths, for 6 layers
% % % 240_42_1
% Obs_BL = [0.41466262,0.4558792,0.50507392,0.50504344,0.51300126,0.4339717,0.44255944,0.39157402,0.38063424,0.39305738,0.38404038,0.37982398,0.33238948,0.31628588;...
% 0.45297344,0.5056505,0.53985922,0.53571648,0.53454554,0.50506122,0.33509204,0.32931354,0.37099494,0.38330632,0.40413178,0.3854831,0.2870835,0.26061416;...
% 0.42030904,0.48281082,0.54699662,0.55175658,0.504698,0.48899572,0.37711126,0.22698456,0.2406523,0.32714438,0.3284601,0.33616646,0.2758186,0.25061672;...
% 0.46190408,0.4895596,0.54720236,0.59442096,0.59171586,0.51962304,0.34720276,0.20550378,0.226695,0.26088594,0.3069463,0.33165542,0.27803094,0.2781173;...
% 0.426847,0.47982378,0.50027586,0.61251846,0.5746496,0.56677052,0.40109902,0.26175208,0.18744438,0.3071368,0.27439874,0.29134562,0.2334895,0.26212292;...
% 0.4185412,0.4694047,0.48911002,0.5334762,0.49486566,0.4750435,0.34272982,0.26179272,0.23569676,0.28867354,0.27715718,0.26702258,0.27318462,0.28847796];
% 
% Mean_MeasuredBL_1 = mean(Obs_BL(1,:))
% Mean_MeasuredBL_2 = mean(Obs_BL(2,:))

%__________________________________________________________________________
% Load Gaussian Process model for K1 and Work done
% GPStruc_Exp_BLmodel_discp = load('gprMdl_Exp');
GPStruc_Model_BLmodel_discp = load('gprMdl_Model_sameStarting');
GP = {GPStruc_Model_BLmodel_discp};

% GPR_Exp = GP{1}.('gprMdl');
GPR_Model = GP{1}.('gprMdl');


% % gp model input requirement: 
% %               T_noz: Nozzle temperature in C, 
% %               v_print: printer speed in mm/sec
% %               x and y coordinates
% T_noz = 240;
% v_print = 42;
% 
% % x & y coordinates of bond lengths
% % Prediction inputs are x&y coordinates of vertical bond (ly) locations
% % x,y coordinate of layer 1 & interface 1 (vertical)
% x = wth*1e3; y = 0.5*hth*1e3; % in mm
% iki_y=y*2;
% TotalLayer=2; % only look at 2 layers, 1st experimental temp data +
% % bond length calculation, 2nd one the heat transfer model + BL model
% Coordinates = zeros(TotalLayer*(numFilaments-1),2);
% for jj=1:TotalLayer
%     for ii=1:numFilaments-1
%         % use x & y coordinates of vertical bonds as predict variables
%         Coordinates((jj-1)*(numFilaments-1) + ii,:) = [ii*x y+(jj-1)*iki_y];
%     end
% end
% 
%     
% Gp_input = [T_noz*ones(TotalLayer*(numFilaments-1),1), ...
%             v_print*ones(TotalLayer*(numFilaments-1),1), ...
%             Coordinates];
% %__________________________________________________________________________
% 
% % % predict BL model discrepency using GP inputs for experimental Temp. data
% % [delta_BL_exp] = predict(GPR_Exp, Gp_input);
% 
% % predict BL model discrepency using GP inputs for model Temp. data
% [delta_BL_model] = predict(GPR_Model, Gp_input);
    



% Material parameters for the bond length model
Surf_tens = 0.29; % Surface tension
b1_ = 0.00345; % Model parameter for temp. dependent surface tension
Visco = 5100; % viscosity
b2_ = 0.0193; % for temp above 130 celcius
% b2_=0.0016;
inputs_Exp = [Surf_tens, b1_, Visco, b2_];


Ref_Temp = 130; %celcius

% inputs_Model = [Surf_tens, b1_, Visco, b2_,alpha];
inputs_Model = [Surf_tens, b1_, Visco];


Kel = 273.15; % 0 Celcius in K
% convert Celcius to Kelvin


% x & y coordinates of bond lengths
% Prediction inputs are x&y coordinates of vertical bond (ly) locations
% x,y coordinate of layer 1 & interface 1 (vertical)
x = wth*1e3; y = 0.5*hth*1e3; % in mm
iki_y=y*2;
% TotalLayer=2; % only look at 2 layers, 1st experimental temp data +
% bond length calculation, 2nd one the heat transfer model + BL model
Coordinates = zeros(1*(numFilaments-1),2);
for jj=1:1 % always obtained model discrepency only for one layer
    for ii=1:numFilaments-1
        % use x & y coordinates of vertical bonds as predict variables
        Coordinates((jj-1)*(numFilaments-1) + ii,:) = [ii*x y+(jj-1)*iki_y];
    end
end




% the layer number which the optimization is carrier out
layerID = 6;
% 'layerID' layers & 15 filaments
Matrix = ones(layerID,15);

% nozzle temperature and printer speed for previous layers that are already
% designed/optimized
% TV: temperature and velocity design variables
% Celcius & m/sec
TV=[259.91 0.043; % 1st layer -> Fval = F1+F2 = -0.3584
    259.65 0.0162; % 2nd layer -> Fval = F1+F2 = -0.3317
    258.63 0.015; % 3rd layer -> Fval = F1+F2 = -0.3176
    258.87 0.015; % 4th layer -> Fval = F1+F2 = -0.3179
    258.99 0.015; % 5th layer -> Fval = F1+F2 = -0.3179
    258.99 0.015]; % 6th layer -> Fval = F1+F2 = -0.3179


% the nozzle temperature is observed to be at least 10 degrees less
% T_N = TV(1)-10; v_p = TV(2);





% Number of samples for Monte Carlo Simulation
n_samples = 30;

time1=tic;
% Run in parallel
poolobj = parpool;

[Fsol,Fval,history] = ProcessDesign_AM_sameStarting_objMeanDominant(Matrix,x_cut,layerID,ProcessParam,...
            inputs_Model,Ref_Temp,TV,GPR_Model,Coordinates, n_samples);

%     delete(poolobj)
toc(time1)







% Sort design variable, printed speed
[sorted_T_N, indices] = sort(history.X(:,1));

% sorted printer speed
sorted_v_p = history.X(indices,2);

% get your output with
sorted_BLs = history.Fval(indices);

% Sorted design variable  
% add 10 to temperature since we already subtracted before to match
% experimental data
sorted_DV = [sorted_T_N+10 sorted_v_p sorted_BLs];
    
    
    
    


% save 30SamplesC_LocalOpt_Miss4
% fid = fopen('Results_30SamplesC_LocalOpt_Miss4.txt', 'wt+');
% format short
% cl = size(history.X, 2);
% for i = 1:size(history.X, 1)
%   fprintf(fid, ' %.3f \t',history.X(i,1:cl));
%   fprintf(fid, ' %.3f \t',history.Fval(i,1));
%   fprintf(fid, ' \n');
% end
% fclose(fid);




