% clc;clear; close all;
% 

% Plot L2 norm
% numTestPoints=floor(0.05*(2:NumTest)*layer*(num_Filaments-1));

% normalize by dividing with number of test points
L2norm = [3.1276e-4,7.8439e-4,0.0017,0.00172845,0.0011855102,.0023381677,.00311329,...
    .003081604,.0018471323, 0.00212521122262, 0.003067793218044, 0.001687814737991,...
    0.002137076110883, 0.002697378605953, 0.002435790005185, 0.001939674095917,...
    0.002956403453499, 0.006686552477393, 0.005255101293957, 0.002797846864382];%./numTestPoints;
% Plot L2 norm of GP model
figure()
loglog(2:NumTest, L2norm, '--');
xlabel('number of samples', 'FontName', 'Times New Roman', ...
       'FontSize',16,'Color','k', 'Interpreter', 'tex');
ylabel('L_2 norm','FontName', 'Times New Roman', ...
       'FontSize',16,'Color','k', 'Interpreter', 'tex');

axis tight;
set(gcf, 'Position',  [200, 300, 900, 600])
set(gcf, 'PaperPosition', [0 0 7 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [7 5]); %Set the paper to have width 5 and height 5.