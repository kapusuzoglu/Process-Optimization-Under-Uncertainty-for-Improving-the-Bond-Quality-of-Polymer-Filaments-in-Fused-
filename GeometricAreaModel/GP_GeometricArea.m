clc; clear; close all;

% Get the Current Directory
Directory = pwd;
% change to current directory
cd(Directory);

% filament dimensions; width & height
rx = 0.8; % width, mm
ry = 0.7; % height, mm

% Average horizontal (inter-layer) bond length
% lx = 0.43; % mm

%%
% input -> GP -> output
% INPUT: TGA (total gap area) 
% OUTPUT: L_b (bond length)
%
% Get the total gap are from the filter algorithm for each sample
TGA = [1.806642, 0.550708576, 0.322521214, 1.077546232, ...
       1.373352, 0.236193]';
L_b = [0.39619555, 0.513391698765432, 0.599268791904762, 0.479333983333333, ...
       0.473316904761905, 0.616576835714286]';
lx = [0.341611, 0.497553, 0.604016, 0.426143];


nl = length(TGA);
X_Lb = TGA; Y_Lb = L_b;






%%
% -------------------------------------------------------------------------
%              Lb/ly (intra-layer bond length)
% -------------------------------------------------------------------------

% randomly select 90% to be training set for GP model
num_train = round(0.9*nl);
idx_Lb = randperm(nl);
% Use the first 90% to train the GP model
idx_Lb_train = idx_Lb(1:num_train);
X_Lb_train = X_Lb(idx_Lb_train, :);
Y_Lb_train = Y_Lb(idx_Lb_train, :);
% The rest will be used to test the GP model
idx_Lb_test = idx_Lb((num_train + 1):end);
X_Lb_test = X_Lb(idx_Lb_test, :);
Y_Lb_test = Y_Lb(idx_Lb_test, :);
% save the training and testing data
save('TrainTestData.mat',...
     'X_Lb_train', 'Y_Lb_train', 'X_Lb_test', 'Y_Lb_test')

% initialize the kernel parameters
sigma0_Lb = std(Y_Lb_train);
sigmaTGA0_Lb = sigma0_Lb;
d = size(X_Lb_train,2);
sigmaM0_Lb = 10*ones(d,1);

% fitting using fitrgp
gprMdl_Lb = fitrgp(X_Lb_train,Y_Lb_train,'Basis','constant',...
                   'FitMethod','exact','PredictMethod','exact',...
                   'KernelFunction','ardsquaredexponential',...
                   'KernelParameters',[sigmaM0_Lb;sigmaTGA0_Lb],...
                   'Sigma',sigma0_Lb,'Standardize',1);
% save the GP model
save('gprMdl_Lb.mat', 'gprMdl_Lb')


% Model parameters
% Trend Function
Beta_Lb= gprMdl_Lb.Beta;
% The kernel parameters. 
sigmaM_Lb = gprMdl_Lb.KernelInformation.KernelParameters(1:end-1,1);
sigmaTGA_Lb = gprMdl_Lb.KernelInformation.KernelParameters(end);
% Estimated noise standard deviation (sigma_n)
sigma_Lb  = gprMdl_Lb.Sigma; 
% save the model parameters
save('Param_gprMdl_Lb.mat',...
     'Beta_Lb', 'sigmaM_Lb', 'sigmaTGA_Lb', 'sigma_Lb')

% test the GP model
Lb_pred = predict(gprMdl_Lb,X_Lb_test);



% -----------------------------
%   Plot GP model test result
% -----------------------------
figure()

plot(1:length(Lb_pred), Y_Lb_test, 'o', 1:length(Lb_pred), Lb_pred, 'x')
legend('True Value', 'GP Prediction')
xlabel('Test point',...
       'interpreter','latex','fontsize', 14)
ylabel('Bond length, $L_b (mm)$',...
       'interpreter','latex','fontsize', 14)
title('Test result - GP Model for $L_b$',...
       'interpreter','latex','fontsize', 16)




% % load array
% loadarray = 1000:1000:8000;
% nl = length(loadarray);
% % crack length array
% clarray = 5:5:90;
% nc = length(clarray);

% % arrange all the data
% X_K1 = zeros(nl*nc, 2); 
% % 1st col: load in lbs, 2nd col: crack length in mm
% Y_K1 = zeros(nl*nc, 1);
% % col: K1 in ( Pa * sqrt(m) )
% 
% k = 1;
% for i = 1:nl % load = loadarray(i) in lbs
%     for j = 1:nc % crack length = clarray(j) in mm
%         filename = ['F_' num2str(loadarray(i)) '_CL_' ...
%                     num2str(clarray(j)) '_K1.txt'];
%         temp = dlmread(filename);
%         Y_K1(k) = temp(1, 3);
%         X_K1(k, 1) = loadarray(i);
%         X_K1(k, 2) = clarray(j);
%         k = k + 1;
%     end
% end


% % randomly select 90% to be training set for GP model
% num_train = round(0.9*nl*nc);
% idx_K1 = randperm(nl*nc);
% % Use the first 90% to train the GP model
% idx_K1_train = idx_K1(1:num_train);
% X_K1_train = X_K1(idx_K1_train, :);
% Y_K1_train = Y_K1(idx_K1_train, :);
% % The rest will be used to test the GP model
% idx_K1_test = idx_K1((num_train + 1):end);
% X_K1_test = X_K1(idx_K1_test, :);
% Y_K1_test = Y_K1(idx_K1_test, :);
% % save the training and testing data
% save('TrainTestData.mat',...
%      'X_K1_train', 'Y_K1_train', 'X_K1_test', 'Y_K1_test')
% 
% % initialize the kernel parameters
% sigma0_K1 = std(Y_K1_train);
% sigmaF0_K1 = sigma0_K1;
% d = size(X_K1_train,2);
% sigmaM0_K1 = 10*ones(d,1);
% 
% % fitting using fitrgp
% gprMdl_K1 = fitrgp(X_K1_train,Y_K1_train,'Basis','constant',...
%                    'FitMethod','exact','PredictMethod','exact',...
%                    'KernelFunction','ardsquaredexponential',...
%                    'KernelParameters',[sigmaM0_K1;sigmaF0_K1],...
%                    'Sigma',sigma0_K1,'Standardize',1);
% % save the GP model
% save('gprMdl_K1.mat', 'gprMdl_K1')
% 
% % Model parameters
% % Trend Function
% Beta_K1= gprMdl_K1.Beta;
% % The kernel parameters. 
% sigmaM_K1 = gprMdl_K1.KernelInformation.KernelParameters(1:end-1,1);
% sigmaF_K1 = gprMdl_K1.KernelInformation.KernelParameters(end);
% % Estimated noise standard deviation (sigma_n)
% sigma_K1  = gprMdl_K1.Sigma; 
% % save the model parameters
% save('Param_gprMdl_K1.mat',...
%      'Beta_K1', 'sigmaM_K1', 'sigmaF_K1', 'sigma_K1')
% 
% % test the GP model
% K1_pred = predict(gprMdl_K1,X_K1_test);
% 
% %%
% % % -----------------------------
% % %   Plot GP model test result
% % % -----------------------------
% % figure()
% % 
% % plot(1:length(K1_pred), Y_K1_test, 'o', 1:length(K1_pred), K1_pred, 'x')
% % legend('True Value', 'GP Prediction')
% % xlabel('Test Point')
% % ylabel('$Stress\ Intensity\ Factor\ K_1\ (MPa \sqrt{m}$)',...
% %        'interpreter','latex','fontsize', 14)
% % title('Test Result - GP Model for K_1','fontsize', 16)
