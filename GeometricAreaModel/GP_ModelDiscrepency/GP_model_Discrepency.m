
tic
%%
% -------------------------------------------------------------------------
%                               GP Model
% -------------------------------------------------------------------------
nl = length(Inp);

% randomly select 95% to be training set for GP model
num_train = round(0.95*nl); 
idx_ = randperm(nl);

% Use the first 95% to train the GP model
idx_train = idx_(1:num_train);
allData_Train = allData(idx_train, :);

% The rest will be used to test the GP model
idx_test = idx_((num_train + 1):end);
allData_Test = allData(idx_test, :);


% % Select test set for GP model and keep it same
% NumGapsIncrement = 14*4;
% GapId = 7; % 7 and 8 are the middle gaps
% Test_set = GapId:NumGapsIncrement:nl;
% all_set = 1:nl;
% rem_set = setdiff(all_set,Test_set);

% 
% % Use these to train the GP model
% % randomly sample some percent of available training data
% Ratio = 0.99;
% Train_set = randsample(rem_set,round(Ratio*length(rem_set)));
% allData_Train = allData(Train_set, :);
% % Test data
% allData_Test = allData(Test_set, :);



% save 'allData_Train.data' allData_Train -ascii
% tbl = readtable('allData_Train.data','Filetype','text',...
% 'ReadVariableNames',false);
% tbl.Properties.VariableNames = {'Temperature','Speed','x','y','ModelDiscrepency'};

% gprMdl = fitrgp(tbl,'ModelDiscrepency','KernelFunction','ardsquaredexponential',...
%       'FitMethod','exact','PredictMethod','exact','Standardize',1);
  
% initialize the kernel parameters
sigma0_Lb = std(allData_Train(:,end));
sigmaTGA0_Lb = sigma0_Lb;
d = size(allData_Train(:,1:end-1),2);
sigmaM0_Lb = 1*ones(d,1);


% fitting using fitrgp
gprMdl = fitrgp(allData_Train(:,1:end-1),allData_Train(:,end),'Basis','constant',...
                   'FitMethod','exact','PredictMethod','exact',...
                   'KernelFunction','ardsquaredexponential',...
                   'KernelParameters',[sigmaM0_Lb;sigmaTGA0_Lb],...
                   'Sigma',sigma0_Lb,'Standardize',1);
               
               
% % fitting using fitrgp
% gprMdl = fitrgp(tbl,'ModelDiscrepency','Basis','constant',...
%                    'FitMethod','exact','PredictMethod','exact',...
%                    'KernelFunction','ardsquaredexponential',...
%                    'KernelParameters',[sigmaM0_Lb;sigmaTGA0_Lb],...
%                    'Sigma',sigma0_Lb,'Standardize',1);

% % Model parameters
% % Trend Function
% Beta_Lb= gprMdl.Beta;
% % The kernel parameters.
% sigmaM_Lb = gprMdl.KernelInformation.KernelParameters(1:end-1,1);
% sigmaTGA_Lb = gprMdl.KernelInformation.KernelParameters(end);
% % Estimated noise standard deviation (sigma_n)
% sigma_Lb  = gprMdl.Sigma;


  
xtest=allData_Test(:,1:end-1);
ytest=allData_Test(:,end);

% test the GP model
[delta_pred, delta_pred_sd, delta_pred_int] = predict(gprMdl,xtest);




% Compute the regression error on the test data
L = loss(gprMdl, xtest, ytest)

d = size(xtest,2);


mean(delta_pred)
std(delta_pred)


% -----------------------------
%   Plot GP model test result
% -----------------------------
figure()
% plot(1:length(delta_pred), ytest, 'or', 1:length(delta_pred), delta_pred, '*b');
% errorbar(1:length(delta_pred), delta_pred,delta_pred_sd,'-b')
boundedline(1:length(delta_pred),delta_pred,1.96*delta_pred_sd, '-bo', 'alpha');hold on;
plot(1:length(delta_pred), ytest, '.m','MarkerSize', 20);
legend('95% Confidence interval','GP predictions', 'Data')
xlabel('Test points', 'FontName', 'Times New Roman', ...
       'FontSize',16,'Color','k', 'Interpreter', 'tex');
ylabel('\delta','FontName', 'Times New Roman', ...
       'FontSize',16,'Color','k', 'Interpreter', 'tex');

axis tight;
set(gcf, 'Position',  [200, 300, 900, 600])
set(gcf, 'PaperPosition', [0 0 7 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [7 5]); %Set the paper to have width 5 and height 5.


toc