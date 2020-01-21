load('gprMdl_K1.mat')
load('Param_gprMdl_K1.mat')
load('TrainTestData.mat')

[K1_train_pred, K1_train_pred_sd] = predict(gprMdl_K1, X_K1_train);
K1_train_pred_ub = K1_train_pred + 1.96 * K1_train_pred_sd;
K1_train_pred_lb = K1_train_pred - 1.96 * K1_train_pred_sd;

[K1_pred, K1_pred_sd] = predict(gprMdl_K1, X_K1_test);
K1_pred_ub = K1_pred + 1.96 * K1_pred_sd;
K1_pred_lb = K1_pred - 1.96 * K1_pred_sd;

figure()
Plot95CBWithGPPrediction(K1_train_pred, K1_train_pred_sd, Y_K1_train)
xlabel('Training Points')
ylabel('$Stress\ Intensity\ Factor\ K_1\ (Pa \sqrt{m}$)',...
       'interpreter','latex','fontsize', 14)
title('Test GP Model for K_1 using Training Points','fontsize', 16)
tightfig;

figure()
Plot95CBWithGPPrediction(K1_pred, K1_pred_sd, Y_K1_test)
xlabel('Testing Points')
ylabel('$Stress\ Intensity\ Factor\ K_1\ (Pa \sqrt{m}$)',...
       'interpreter','latex','fontsize', 14)
title('Test GP Model for K_1 using Testing Points','fontsize', 16)
tightfig;

function Plot95CBWithGPPrediction(Data_GPPred, SD, Data)
% 
% Input Requirement
% Data_GPPred: Array of mean prediction from GP Model
% SD: standard deviation associated with the GP Prediction
% Data: Real value of the data
% 
% Total number of data
n = length(Data);
% Upper bound for 95 confidence interval
UB = Data_GPPred + 1.96 * SD;
% Lower
LB = Data_GPPred - 1.96 * SD;

% Plot actual value of the data
plot(1:n, Data, 'xr'); hold on
% Plot the GP prediction
plot(1:n, Data_GPPred, 'ok'); hold on;
% Plot 95% confidence bound
for i = 1:n
    plot([i, i], [UB(i), LB(i)], 'k'); hold on
end

hold off
legend('True Value (Abaqus)', 'GP Model Prediction with 95% Confidence Bound')
% xlabel('$Training\ Points$', 'interpreter','latex','fontsize', 14)
% ylabel('$Stress\ Intensity\ Factor\ K_1\ (MPa \sqrt{m}$)',...
%        'interpreter','latex','fontsize', 14)
% title('Test GP Model for K_1 using Training Points','fontsize', 16)
% set(gca,'xtick',[])
end