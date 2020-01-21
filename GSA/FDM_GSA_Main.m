%% Load the test data
clear; clc;
close all;
% rng(1)


% ======================================================================= %
% ======================================================================= %
% ===========================      GSA     ============================== %
% ======================================================================= %
% ======================================================================= %

% n: the number of input samples you'd like to use in each loop. For
% example, N can be 10000 but only 5000 are used in the analysis if you
% set n=5000
% N = 4;
% n = N^3-5*N;


KK = 273.15; % 0 Celcius in K

%Process Variables
%________________________________________________________________________
T_L = 215; %Extrusion temperature (ºC)
T_E = 110; %Temperature of the envelope (ºC)
v = 0.032; %Velocity of the extrusion head (m/sec)
w = 0.8e-3; % layer width [meters]
h = 0.7e-3; % layer height [meters]
L = 0.035; %Length of the filament (meters)
numLayers = 6;
number_filament = 15; 

% Process parameters
ProcessParam = [T_L, T_E, v, w, h, L];



% =======================         MC Samples      ======================= % 
N_Samples = 50;
% icdf of the parameters
U = lhsdesign(N_Samples, 4);

% Material Properties
%________________________________________________________________________
%Density (kg/m^3)
ro_mu = 1040; ro_cov = 0.05; % material A
ro_Samples = icdf('Normal', U(:,1), ro_mu, ro_cov*ro_mu);
%Specific heat (J/kg.K)
C_mu = 1290; C_cov = 0.05; % material A
C_Samples = icdf('Normal', U(:,2), C_mu, C_cov*C_mu);
% Heat transfer coefficient (loss of heat by natural convection)
h_conv_mu = 86; h_conv_cov = 0.05;
hconv_Samples = icdf('Normal', U(:,3), h_conv_mu, h_conv_cov*h_conv_mu);

% % %Thermal contact conductances between
% h_cond_mu = 86; h_cond_cov = 0.1; % filament and adjacent filament & support
% hcond_Samples = icdf('Normal', U(:,4), h_cond_mu, h_cond_cov*h_cond_mu);


% %Thermal contact conductances between
lam_mu = .2; lam_cov = 0.1; % surface contact fractions
muN_lam = log(lam_mu/sqrt(1 + lam_cov^2));
sigmaN_lam = sqrt(log(1 + lam_cov^2));
lam_Samples = icdf('LogNormal', U(:,4), muN_lam, sigmaN_lam);



% Material parameters
MaterialParam = [ro_Samples, C_Samples, hconv_Samples, lam_Samples];



% n = 3;
% tic
% %% Double loop MC method
% First_order_double_loop   = DoubleLoopGSA_First(input_sample1, @Forman_Model_5000, n, Other_inputs);
% Total_effects_double_loop = DoubleLoopGSA_Total(input_sample1, @Forman_Model_5000, n, Other_inputs);
% toc




tic
% ======================================================================= %
%% Costa Model - Analytical
Matrix = ones(numLayers,number_filament);
x = L/2; % location of the cut in meters

% [Temp_C,time,limite] = FDM_TempModel_GSA(Matrix, x, ProcessParam, MaterialParam);

Temp_C2=zeros(N_Samples,1006);
Temp_C3=zeros(N_Samples,1006);
Temp_C4=zeros(N_Samples,1006);
Temp_C87=zeros(N_Samples,1006);
Temp_C88=zeros(N_Samples,1006);
Temp_C89=zeros(N_Samples,1006);

ticBytes(gcp);
parfor ii=1:N_Samples
    [Temp_C2(ii,:),Temp_C3(ii,:),Temp_C4(ii,:), Temp_C87(ii,:),Temp_C88(ii,:),Temp_C89(ii,:)] = FDM_TempModel_GSA(Matrix, x, ProcessParam, MaterialParam(ii,:));
end
Temp_C2=Temp_C2(:,17:end);
Temp_C3=Temp_C3(:,27:end);
Temp_C4=Temp_C4(:,37:end);

Temp_C87=Temp_C87(:,867:end);
Temp_C88=Temp_C88(:,877:end);
Temp_C89=Temp_C89(:,887:end);

tocBytes(gcp)

% % convert Celcius to Kelvin
% Temp = Temp_C + KK;
% 
% % Average two adjacent lines' temperatures to obtain interface temperatures
% % for each interface on that layer
% int_Temp = zeros(size(Temp,1),numFilaments-1);
% T_idx = zeros(1,numFilaments-1);
% for ii=2:numFilaments
%     T_idx(ii-1) = find(Temp(:,ii)<T_L+KK, 1, 'first');
%     int_Temp(T_idx(ii-1):end,ii-1) = (Temp(T_idx(ii-1):end,ii-1)+Temp(T_idx(ii-1):end,ii))/2;
% end
% 
% t_final=time(end);
% 
% 
% % Neck growth calculation using the temperature field from Costa model
% for rr=1:numFilaments-1 % loop over each interface for layer 1
%     num=size(int_Temp,1)-T_idx(rr)+1;
%     dt=t_final/num;
%     [bond3,t_bond3,theta3]=Neck_Growth_Model(int_Temp(T_idx(rr):end,rr),dt,t_final,w,h);
%     % method 2, more advisable
%     my_bonds = strcat('bond3_',num2str(rr));
%     variable.(my_bonds) = bond3;
%     my_bonds = strcat('t_bond3_',num2str(rr));
%     variable.(my_bonds) = t_bond3;
%     my_bonds = strcat('theta3_',num2str(rr));
%     variable.(my_bonds) = theta3;
% end

% ======================================================================= %
fprintf('\nAnalytical model is over: ');toc
fprintf('\n');


tic
%% Single loop methods
% CDF

input_cdf = {[1 2 3 4], {'Normal','Normal','Normal','LogNormal'}, ...
    {ro_mu,ro_cov*ro_mu}, {C_mu,C_cov*C_mu}, {h_conv_mu,h_conv_cov*h_conv_mu}, {muN_lam, sigmaN_lam}};

SingleFirst_Index_Alg1_fil2 = SingleLoop_FirstOrder_Alg1(MaterialParam, ...
    Temp_C2, ceil(sqrt(N_Samples)),  input_cdf);
SingleFirst_Index_Alg2_fil2 = SingleLoop_FirstOrder_Alg2(MaterialParam, ...
    Temp_C2, ceil(sqrt(N_Samples)),  input_cdf);

fprintf('Sum of indices: %f\n', sum(SingleFirst_Index_Alg1_fil2(1,:)));
fprintf('Sum of indices: %f \n', sum(SingleFirst_Index_Alg2_fil2(1,:)));

%% Single loop methods
SingleFirst_Index_Alg1_fil3 = SingleLoop_FirstOrder_Alg1(MaterialParam, ...
    Temp_C3, ceil(sqrt(N_Samples)),  input_cdf);
SingleFirst_Index_Alg2_fil3 = SingleLoop_FirstOrder_Alg2(MaterialParam, ...
    Temp_C3, ceil(sqrt(N_Samples)),  input_cdf);
fprintf('Sum of indices: %f\n', sum(SingleFirst_Index_Alg1_fil3(1,:)));
fprintf('Sum of indices: %f \n', sum(SingleFirst_Index_Alg2_fil3(1,:)));

%% Single loop methods
SingleFirst_Index_Alg1_fil4 = SingleLoop_FirstOrder_Alg1(MaterialParam, ...
    Temp_C4, ceil(sqrt(N_Samples)),  input_cdf);
SingleFirst_Index_Alg2_fil4 = SingleLoop_FirstOrder_Alg2(MaterialParam, ...
    Temp_C4, ceil(sqrt(N_Samples)),  input_cdf);
fprintf('Sum of indices: %f\n', sum(SingleFirst_Index_Alg1_fil4(1,:)));
fprintf('Sum of indices: %f\n', sum(SingleFirst_Index_Alg2_fil4(1,:)));

%% Single loop methods
SingleFirst_Index_Alg1_fil87 = SingleLoop_FirstOrder_Alg1(MaterialParam, ...
    Temp_C87, ceil(sqrt(N_Samples)),  input_cdf);
SingleFirst_Index_Alg2_fil87 = SingleLoop_FirstOrder_Alg2(MaterialParam, ...
    Temp_C87, ceil(sqrt(N_Samples)),  input_cdf);
fprintf('Sum of indices: %f\n', sum(SingleFirst_Index_Alg1_fil87(1,:)));
fprintf('Sum of indices: %f\n', sum(SingleFirst_Index_Alg2_fil87(1,:)));

%% Single loop methods
SingleFirst_Index_Alg1_fil88 = SingleLoop_FirstOrder_Alg1(MaterialParam, ...
    Temp_C88, ceil(sqrt(N_Samples)),  input_cdf);
SingleFirst_Index_Alg2_fil88 = SingleLoop_FirstOrder_Alg2(MaterialParam, ...
    Temp_C88, ceil(sqrt(N_Samples)),  input_cdf);

fprintf('Sum of indices: %f\n', sum(SingleFirst_Index_Alg1_fil88(1,:)));
fprintf('Sum of indices: %f\n', sum(SingleFirst_Index_Alg2_fil88(1,:)));

%% Single loop methods

SingleFirst_Index_Alg1_fil89 = SingleLoop_FirstOrder_Alg1(MaterialParam, ...
    Temp_C89, ceil(sqrt(N_Samples)),  input_cdf);
SingleFirst_Index_Alg2_fil89 = SingleLoop_FirstOrder_Alg2(MaterialParam, ...
    Temp_C89, ceil(sqrt(N_Samples)),  input_cdf);

fprintf('Sum of indices: %f\n', sum(SingleFirst_Index_Alg1_fil89(1,:)));
fprintf('Sum of indices: %f\n', sum(SingleFirst_Index_Alg2_fil89(1,:)));



fprintf('==========     Single Loop completed    ========\n');

toc



% 
% N = 10;
% %% Latin Hypercube Samples
% Mu = [ro_mu,C_mu,h_conv_mu,h_cond_mu]; % mean value
% Std = Mu.*[ro_cov C_cov h_conv_cov h_cond_cov]; % standard deviation
% % Sigma = diag(Std.^2); % standard deviation;
% 
% % input_sample = lhsnorm(Mu, Std, N);
% % A = input_sample;
% % B = lhsnorm(Mu, Std, N);
% 
% k = size(Mu, 2);
% Other_inputs = [x, ProcessParam];
% other_inputs = [Mu, Std, Other_inputs];
% First_order_FAST = Sen_FirstOrder_RBD(N, k, @Temp_Model_RBD,other_inputs);
% fprintf('\n Improved FAST completed: '); toc




t=1:size(SingleFirst_Index_Alg1_fil2,1);
t=t/10;

% plot

% fil. 2
% fname = 'C:\Users\kapusub\Documents\MATLAB\TermPaper\GSA_Prognosis032';
figure
set(0,'DefaultLineLineWidth',2)
plot(t,SingleFirst_Index_Alg1_fil2(:,1), 'm-'); hold on;
plot(t,SingleFirst_Index_Alg2_fil2(:,1), 'b--'); hold on;
ylim([0, 1]);
xlabel('Time (s)', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
ylabel('First-order index', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
legend({'Single loop Algorithm 1', 'Single loop Algorithm 2'},'location', 'northwest', 'FontName', 'Times New Roman', ...
'FontSize',10,'Interpreter', 'LaTeX')
title('$\rho$','FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
set(gcf, 'Position',  [200, 300, 900, 600])
set(gcf, 'PaperPosition', [0 0 7 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [7 5]); %Set the paper to have width 5 and height 5.
% saveas(gcf, fullfile(fname,'\Fig_UniMaxF_3_1'), 'pdf') %Save figure


figure
set(0,'DefaultLineLineWidth',2)
plot(t,SingleFirst_Index_Alg1_fil2(:,2), 'm-'); hold on;
plot(t,SingleFirst_Index_Alg2_fil2(:,2), 'b--'); hold on;
ylim([0, 1]);
xlabel('Time (s)', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
ylabel('First-order index', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
legend({'Single loop Algorithm 1', 'Single loop Algorithm 2'},'location', 'northwest', 'FontName', 'Times New Roman', ...
'FontSize',10,'Interpreter', 'LaTeX')
title('C','FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
set(gcf, 'Position',  [200, 300, 900, 600])
set(gcf, 'PaperPosition', [0 0 7 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [7 5]); %Set the paper to have width 5 and height 5.
% saveas(gcf, fullfile(fname,'\Fig_UniMaxF_3_1'), 'pdf') %Save figure


figure
set(0,'DefaultLineLineWidth',2)
plot(t,SingleFirst_Index_Alg1_fil2(:,3), 'm-'); hold on;
plot(t,SingleFirst_Index_Alg2_fil2(:,3), 'b--'); hold on;
ylim([0, 1]);
xlabel('Time (s)', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
ylabel('First-order index', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
legend({'Single loop Algorithm 1', 'Single loop Algorithm 2'},'location', 'northwest', 'FontName', 'Times New Roman', ...
'FontSize',10,'Interpreter', 'LaTeX')
title('$h_{conv}$','FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
set(gcf, 'Position',  [200, 300, 900, 600])
set(gcf, 'PaperPosition', [0 0 7 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [7 5]); %Set the paper to have width 5 and height 5.
% saveas(gcf, fullfile(fname,'\Fig_UniMaxF_3_1'), 'pdf') %Save figure

figure
set(0,'DefaultLineLineWidth',2)
plot(t,SingleFirst_Index_Alg1_fil2(:,4), 'm-'); hold on;
plot(t,SingleFirst_Index_Alg2_fil2(:,4), 'b--'); hold on;
ylim([0, 1]);
xlabel('Time (s)', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
ylabel('First-order index', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
legend({'Single loop Algorithm 1', 'Single loop Algorithm 2'},'location', 'northwest', 'FontName', 'Times New Roman', ...
'FontSize',10,'Interpreter', 'LaTeX')
title('$h_{cond}$','FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
set(gcf, 'Position',  [200, 300, 900, 600])
set(gcf, 'PaperPosition', [0 0 7 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [7 5]); %Set the paper to have width 5 and height 5.
% saveas(gcf, fullfile(fname,'\Fig_UniMaxF_3_1'), 'pdf') %Save figure

t=1:size(SingleFirst_Index_Alg1_fil3,1);
t=t/10;
% fil. 3
figure
set(0,'DefaultLineLineWidth',2)
plot(t,SingleFirst_Index_Alg1_fil3(:,1), 'm-'); hold on;
plot(t,SingleFirst_Index_Alg2_fil3(:,1), 'b--'); hold on;
ylim([0, 1]);
xlabel('Time (s)', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
ylabel('First-order index', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
legend({'Single loop Algorithm 1', 'Single loop Algorithm 2'},'location', 'northwest', 'FontName', 'Times New Roman', ...
'FontSize',10,'Interpreter', 'LaTeX')
title('$\rho$','FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
set(gcf, 'Position',  [200, 300, 900, 600])
set(gcf, 'PaperPosition', [0 0 7 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [7 5]); %Set the paper to have width 5 and height 5.
% saveas(gcf, fullfile(fname,'\Fig_UniMaxF_3_1'), 'pdf') %Save figure


figure
set(0,'DefaultLineLineWidth',2)
plot(t,SingleFirst_Index_Alg1_fil3(:,2), 'm-'); hold on;
plot(t,SingleFirst_Index_Alg2_fil3(:,2), 'b--'); hold on;
ylim([0, 1]);
xlabel('Time (s)', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
ylabel('First-order index', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
legend({'Single loop Algorithm 1', 'Single loop Algorithm 2'},'location', 'northwest', 'FontName', 'Times New Roman', ...
'FontSize',10,'Interpreter', 'LaTeX')
title('C','FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
set(gcf, 'Position',  [200, 300, 900, 600])
set(gcf, 'PaperPosition', [0 0 7 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [7 5]); %Set the paper to have width 5 and height 5.
% saveas(gcf, fullfile(fname,'\Fig_UniMaxF_3_1'), 'pdf') %Save figure


figure
set(0,'DefaultLineLineWidth',2)
plot(t,SingleFirst_Index_Alg1_fil3(:,3), 'm-'); hold on;
plot(t,SingleFirst_Index_Alg2_fil3(:,3), 'b--'); hold on;
ylim([0, 1]);
xlabel('Time (s)', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
ylabel('First-order index', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
legend({'Single loop Algorithm 1', 'Single loop Algorithm 2'},'location', 'northwest', 'FontName', 'Times New Roman', ...
'FontSize',10,'Interpreter', 'LaTeX')
title('$h_{conv}$','FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
set(gcf, 'Position',  [200, 300, 900, 600])
set(gcf, 'PaperPosition', [0 0 7 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [7 5]); %Set the paper to have width 5 and height 5.
% saveas(gcf, fullfile(fname,'\Fig_UniMaxF_3_1'), 'pdf') %Save figure

figure
set(0,'DefaultLineLineWidth',2)
plot(t,SingleFirst_Index_Alg1_fil3(:,4), 'm-'); hold on;
plot(t,SingleFirst_Index_Alg2_fil3(:,4), 'b--'); hold on;
ylim([0, 1]);
xlabel('Time (s)', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
ylabel('First-order index', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
legend({'Single loop Algorithm 1', 'Single loop Algorithm 2'},'location', 'northwest', 'FontName', 'Times New Roman', ...
'FontSize',10,'Interpreter', 'LaTeX')
title('$h_{cond}$','FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
set(gcf, 'Position',  [200, 300, 900, 600])
set(gcf, 'PaperPosition', [0 0 7 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [7 5]); %Set the paper to have width 5 and height 5.
% saveas(gcf, fullfile(fname,'\Fig_UniMaxF_3_1'), 'pdf') %Save figure


t=1:size(SingleFirst_Index_Alg1_fil4,1);
t=t/10;
% fil. 4
figure
set(0,'DefaultLineLineWidth',2)
plot(t,SingleFirst_Index_Alg1_fil4(:,1), 'm-'); hold on;
plot(t,SingleFirst_Index_Alg2_fil4(:,1), 'b--'); hold on;
ylim([0, 1]);
xlabel('Time (s)', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
ylabel('First-order index', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
legend({'Single loop Algorithm 1', 'Single loop Algorithm 2'},'location', 'northwest', 'FontName', 'Times New Roman', ...
'FontSize',10,'Interpreter', 'LaTeX')
title('$\rho$','FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
set(gcf, 'Position',  [200, 300, 900, 600])
set(gcf, 'PaperPosition', [0 0 7 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [7 5]); %Set the paper to have width 5 and height 5.
% saveas(gcf, fullfile(fname,'\Fig_UniMaxF_3_1'), 'pdf') %Save figure


figure
set(0,'DefaultLineLineWidth',2)
plot(t,SingleFirst_Index_Alg1_fil4(:,2), 'm-'); hold on;
plot(t,SingleFirst_Index_Alg2_fil4(:,2), 'b--'); hold on;
ylim([0, 1]);
xlabel('Time (s)', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
ylabel('First-order index', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
legend({'Single loop Algorithm 1', 'Single loop Algorithm 2'},'location', 'northwest', 'FontName', 'Times New Roman', ...
'FontSize',10,'Interpreter', 'LaTeX')
title('C','FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
set(gcf, 'Position',  [200, 300, 900, 600])
set(gcf, 'PaperPosition', [0 0 7 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [7 5]); %Set the paper to have width 5 and height 5.
% saveas(gcf, fullfile(fname,'\Fig_UniMaxF_3_1'), 'pdf') %Save figure


figure
set(0,'DefaultLineLineWidth',2)
plot(t,SingleFirst_Index_Alg1_fil4(:,3), 'm-'); hold on;
plot(t,SingleFirst_Index_Alg2_fil4(:,3), 'b--'); hold on;
ylim([0, 1]);
xlabel('Time (s)', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
ylabel('First-order index', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
legend({'Single loop Algorithm 1', 'Single loop Algorithm 2'},'location', 'northwest', 'FontName', 'Times New Roman', ...
'FontSize',10,'Interpreter', 'LaTeX')
title('$h_{conv}$','FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
set(gcf, 'Position',  [200, 300, 900, 600])
set(gcf, 'PaperPosition', [0 0 7 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [7 5]); %Set the paper to have width 5 and height 5.
% saveas(gcf, fullfile(fname,'\Fig_UniMaxF_3_1'), 'pdf') %Save figure

figure
set(0,'DefaultLineLineWidth',2)
plot(t,SingleFirst_Index_Alg1_fil4(:,4), 'm-'); hold on;
plot(t,SingleFirst_Index_Alg2_fil4(:,4), 'b--'); hold on;
ylim([0, 1]);
xlabel('Time (s)', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
ylabel('First-order index', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
legend({'Single loop Algorithm 1', 'Single loop Algorithm 2'},'location', 'northwest', 'FontName', 'Times New Roman', ...
'FontSize',10,'Interpreter', 'LaTeX')
title('$h_{cond}$','FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
set(gcf, 'Position',  [200, 300, 900, 600])
set(gcf, 'PaperPosition', [0 0 7 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [7 5]); %Set the paper to have width 5 and height 5.
% saveas(gcf, fullfile(fname,'\Fig_UniMaxF_3_1'), 'pdf') %Save figure

















t=1:size(SingleFirst_Index_Alg1_fil87,1);
t=t/10;


% fil. 87
figure
set(0,'DefaultLineLineWidth',2)
plot(t,SingleFirst_Index_Alg1_fil87(:,1), 'm-'); hold on;
plot(t,SingleFirst_Index_Alg2_fil87(:,1), 'b--'); hold on;
ylim([0, 1]);
xlabel('Time (s)', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
ylabel('First-order index', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
legend({'Single loop Algorithm 1', 'Single loop Algorithm 2'},'location', 'northwest', 'FontName', 'Times New Roman', ...
'FontSize',10,'Interpreter', 'LaTeX')
title('$\rho$','FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
set(gcf, 'Position',  [200, 300, 900, 600])
set(gcf, 'PaperPosition', [0 0 7 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [7 5]); %Set the paper to have width 5 and height 5.
% saveas(gcf, fullfile(fname,'\Fig_UniMaxF_3_1'), 'pdf') %Save figure


figure
set(0,'DefaultLineLineWidth',2)
plot(t,SingleFirst_Index_Alg1_fil87(:,2), 'm-'); hold on;
plot(t,SingleFirst_Index_Alg2_fil87(:,2), 'b--'); hold on;
ylim([0, 1]);
xlabel('Time (s)', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
ylabel('First-order index', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
legend({'Single loop Algorithm 1', 'Single loop Algorithm 2'},'location', 'northwest', 'FontName', 'Times New Roman', ...
'FontSize',10,'Interpreter', 'LaTeX')
title('C','FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
set(gcf, 'Position',  [200, 300, 900, 600])
set(gcf, 'PaperPosition', [0 0 7 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [7 5]); %Set the paper to have width 5 and height 5.
% saveas(gcf, fullfile(fname,'\Fig_UniMaxF_3_1'), 'pdf') %Save figure


figure
set(0,'DefaultLineLineWidth',2)
plot(t,SingleFirst_Index_Alg1_fil87(:,3), 'm-'); hold on;
plot(t,SingleFirst_Index_Alg2_fil87(:,3), 'b--'); hold on;
ylim([0, 1]);
xlabel('Time (s)', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
ylabel('First-order index', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
legend({'Single loop Algorithm 1', 'Single loop Algorithm 2'},'location', 'northwest', 'FontName', 'Times New Roman', ...
'FontSize',10,'Interpreter', 'LaTeX')
title('$h_{conv}$','FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
set(gcf, 'Position',  [200, 300, 900, 600])
set(gcf, 'PaperPosition', [0 0 7 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [7 5]); %Set the paper to have width 5 and height 5.
% saveas(gcf, fullfile(fname,'\Fig_UniMaxF_3_1'), 'pdf') %Save figure

figure
set(0,'DefaultLineLineWidth',2)
plot(t,SingleFirst_Index_Alg1_fil87(:,4), 'm-'); hold on;
plot(t,SingleFirst_Index_Alg2_fil87(:,4), 'b--'); hold on;
ylim([0, 1]);
xlabel('Time (s)', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
ylabel('First-order index', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
legend({'Single loop Algorithm 1', 'Single loop Algorithm 2'},'location', 'northwest', 'FontName', 'Times New Roman', ...
'FontSize',10,'Interpreter', 'LaTeX')
title('$h_{cond}$','FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
set(gcf, 'Position',  [200, 300, 900, 600])
set(gcf, 'PaperPosition', [0 0 7 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [7 5]); %Set the paper to have width 5 and height 5.
% saveas(gcf, fullfile(fname,'\Fig_UniMaxF_3_1'), 'pdf') %Save figure

t=1:size(SingleFirst_Index_Alg1_fil88,1);
t=t/10;
% fil. 88
figure
set(0,'DefaultLineLineWidth',2)
plot(t,SingleFirst_Index_Alg1_fil88(:,1), 'm-'); hold on;
plot(t,SingleFirst_Index_Alg2_fil88(:,1), 'b--'); hold on;
ylim([0, 1]);
xlabel('Time (s)', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
ylabel('First-order index', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
legend({'Single loop Algorithm 1', 'Single loop Algorithm 2'},'location', 'northwest', 'FontName', 'Times New Roman', ...
'FontSize',10,'Interpreter', 'LaTeX')
title('$\rho$','FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
set(gcf, 'Position',  [200, 300, 900, 600])
set(gcf, 'PaperPosition', [0 0 7 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [7 5]); %Set the paper to have width 5 and height 5.
% saveas(gcf, fullfile(fname,'\Fig_UniMaxF_3_1'), 'pdf') %Save figure


figure
set(0,'DefaultLineLineWidth',2)
plot(t,SingleFirst_Index_Alg1_fil88(:,2), 'm-'); hold on;
plot(t,SingleFirst_Index_Alg2_fil88(:,2), 'b--'); hold on;
ylim([0, 1]);
xlabel('Time (s)', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
ylabel('First-order index', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
legend({'Single loop Algorithm 1', 'Single loop Algorithm 2'},'location', 'northwest', 'FontName', 'Times New Roman', ...
'FontSize',10,'Interpreter', 'LaTeX')
title('C','FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
set(gcf, 'Position',  [200, 300, 900, 600])
set(gcf, 'PaperPosition', [0 0 7 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [7 5]); %Set the paper to have width 5 and height 5.
% saveas(gcf, fullfile(fname,'\Fig_UniMaxF_3_1'), 'pdf') %Save figure


figure
set(0,'DefaultLineLineWidth',2)
plot(t,SingleFirst_Index_Alg1_fil88(:,3), 'm-'); hold on;
plot(t,SingleFirst_Index_Alg2_fil88(:,3), 'b--'); hold on;
ylim([0, 1]);
xlabel('Time (s)', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
ylabel('First-order index', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
legend({'Single loop Algorithm 1', 'Single loop Algorithm 2'},'location', 'northwest', 'FontName', 'Times New Roman', ...
'FontSize',10,'Interpreter', 'LaTeX')
title('$h_{conv}$','FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
set(gcf, 'Position',  [200, 300, 900, 600])
set(gcf, 'PaperPosition', [0 0 7 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [7 5]); %Set the paper to have width 5 and height 5.
% saveas(gcf, fullfile(fname,'\Fig_UniMaxF_3_1'), 'pdf') %Save figure

figure
set(0,'DefaultLineLineWidth',2)
plot(t,SingleFirst_Index_Alg1_fil88(:,4), 'm-'); hold on;
plot(t,SingleFirst_Index_Alg2_fil88(:,4), 'b--'); hold on;
ylim([0, 1]);
xlabel('Time (s)', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
ylabel('First-order index', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
legend({'Single loop Algorithm 1', 'Single loop Algorithm 2'},'location', 'northwest', 'FontName', 'Times New Roman', ...
'FontSize',10,'Interpreter', 'LaTeX')
title('$h_{cond}$','FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
set(gcf, 'Position',  [200, 300, 900, 600])
set(gcf, 'PaperPosition', [0 0 7 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [7 5]); %Set the paper to have width 5 and height 5.
% saveas(gcf, fullfile(fname,'\Fig_UniMaxF_3_1'), 'pdf') %Save figure


t=1:size(SingleFirst_Index_Alg1_fil89,1);
t=t/10;
% fil. 89
figure
set(0,'DefaultLineLineWidth',2)
plot(t,SingleFirst_Index_Alg1_fil89(:,1), 'm-'); hold on;
plot(t,SingleFirst_Index_Alg2_fil89(:,1), 'b--'); hold on;
ylim([0, 1]);
xlabel('Time (s)', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
ylabel('First-order index', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
legend({'Single loop Algorithm 1', 'Single loop Algorithm 2'},'location', 'northwest', 'FontName', 'Times New Roman', ...
'FontSize',10,'Interpreter', 'LaTeX')
title('$\rho$','FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
set(gcf, 'Position',  [200, 300, 900, 600])
set(gcf, 'PaperPosition', [0 0 7 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [7 5]); %Set the paper to have width 5 and height 5.
% saveas(gcf, fullfile(fname,'\Fig_UniMaxF_3_1'), 'pdf') %Save figure


figure
set(0,'DefaultLineLineWidth',2)
plot(t,SingleFirst_Index_Alg1_fil89(:,2), 'm-'); hold on;
plot(t,SingleFirst_Index_Alg2_fil89(:,2), 'b--'); hold on;
ylim([0, 1]);
xlabel('Time (s)', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
ylabel('First-order index', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
legend({'Single loop Algorithm 1', 'Single loop Algorithm 2'},'location', 'northwest', 'FontName', 'Times New Roman', ...
'FontSize',10,'Interpreter', 'LaTeX')
title('C','FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
set(gcf, 'Position',  [200, 300, 900, 600])
set(gcf, 'PaperPosition', [0 0 7 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [7 5]); %Set the paper to have width 5 and height 5.
% saveas(gcf, fullfile(fname,'\Fig_UniMaxF_3_1'), 'pdf') %Save figure


figure
set(0,'DefaultLineLineWidth',2)
plot(t,SingleFirst_Index_Alg1_fil89(:,3), 'm-'); hold on;
plot(t,SingleFirst_Index_Alg2_fil89(:,3), 'b--'); hold on;
ylim([0, 1]);
xlabel('Time (s)', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
ylabel('First-order index', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
legend({'Single loop Algorithm 1', 'Single loop Algorithm 2'},'location', 'northwest', 'FontName', 'Times New Roman', ...
'FontSize',10,'Interpreter', 'LaTeX')
title('$h_{conv}$','FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
set(gcf, 'Position',  [200, 300, 900, 600])
set(gcf, 'PaperPosition', [0 0 7 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [7 5]); %Set the paper to have width 5 and height 5.
% saveas(gcf, fullfile(fname,'\Fig_UniMaxF_3_1'), 'pdf') %Save figure


figure
set(0,'DefaultLineLineWidth',2)
plot(t,SingleFirst_Index_Alg1_fil89(:,4), 'm-'); hold on;
plot(t,SingleFirst_Index_Alg2_fil89(:,4), 'b--'); hold on;
ylim([0, 1]);
xlabel('Time (s)', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
ylabel('First-order index', 'FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
legend({'Single loop Algorithm 1', 'Single loop Algorithm 2'},'location', 'northwest', 'FontName', 'Times New Roman', ...
'FontSize',10,'Interpreter', 'LaTeX')
title('$h_{cond}$','FontName', 'Times New Roman', ...
'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
set(gcf, 'Position',  [200, 300, 900, 600])
set(gcf, 'PaperPosition', [0 0 7 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [7 5]); %Set the paper to have width 5 and height 5.
% saveas(gcf, fullfile(fname,'\Fig_UniMaxF_3_1'), 'pdf') %Save figure






























% 
% %% plot
% fname = 'C:\Users\kapusub\Documents\MATLAB\TermPaper\GSA_Prognosis032';
% figure
% set(0,'DefaultLineLineWidth',2)
% % plot(First_order_double_loop, 'r+');
% % hold on;
% plot(t,SingleFirst_Index_Alg1(:,1), 'm-'); hold on;
% plot(t,SingleFirst_Index_Alg2(:,1), 'b--'); hold on;
% % plot(t,First_order_Saltelli(2:end,1), 'c.-.'); hold on;
% % plot(t,First_order_FAST(:,1), 'r-.'); hold on;
% % plot(t,First_order_Sobol(2:end,1), 'k:');
% ylim([0, 1]);
% xlabel('Time (Number of cycles)', 'FontName', 'Times New Roman', ...
% 'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
% ylabel('First-order index', 'FontName', 'Times New Roman', ...
% 'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
% legend({'Single loop Algorithm 1', 'Single loop Algorithm 2'},'location', 'northwest', 'FontName', 'Times New Roman', ...
% 'FontSize',10,'Interpreter', 'LaTeX')
% title('C','FontName', 'Times New Roman', ...
% 'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
% set(gcf, 'Position',  [200, 300, 900, 600])
% set(gcf, 'PaperPosition', [0 0 7 5]); %Position plot at left hand corner with width 5 and height 5.
% set(gcf, 'PaperSize', [7 5]); %Set the paper to have width 5 and height 5.
% saveas(gcf, fullfile(fname,'\Fig_UniMaxF_3_1'), 'pdf') %Save figure
% %% plot
% figure
% % plot(First_order_double_loop, 'r+');
% % hold on;
% plot(t,SingleFirst_Index_Alg1(:,2), 'm-'); hold on;
% plot(t,SingleFirst_Index_Alg2(:,2), 'b--'); hold on;
% % plot(t,First_order_Saltelli(2:end,2), 'c.-.'); hold on;
% % plot(t,First_order_FAST(:,2), 'r-.'); hold on;
% % plot(t,First_order_Sobol(2:end,2), 'k:');
% ylim([0, 1]);
% xlabel('Time (Number of cycles)', 'FontName', 'Times New Roman', ...
% 'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
% ylabel('First-order index', 'FontName', 'Times New Roman', ...
% 'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
% legend({'Single loop Algorithm 1', 'Single loop Algorithm 2'},'location', 'northwest', 'FontName', 'Times New Roman', ...
% 'FontSize',10,'Interpreter', 'LaTeX')
% title('m','FontName', 'Times New Roman', ...
% 'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
% set(gcf, 'Position',  [200, 300, 900, 600])
% set(gcf, 'PaperPosition', [0 0 7 5]); %Position plot at left hand corner with width 5 and height 5.
% set(gcf, 'PaperSize', [7 5]); %Set the paper to have width 5 and height 5.
% saveas(gcf, fullfile(fname,'\Fig_UniMaxF_3_2'), 'pdf') %Save figure
% %% plot
% figure
% % plot(First_order_double_loop, 'r+');
% % hold on;
% plot(t,SingleFirst_Index_Alg1(:,3), 'm-'); hold on;
% plot(t,SingleFirst_Index_Alg2(:,3), 'b--'); hold on;
% % plot(t,First_order_Saltelli(2:end,3), 'c.-.'); hold on;
% % plot(t,First_order_FAST(:,3), 'r-.'); hold on;
% % plot(t,First_order_Sobol(2:end,3), 'k:');
% ylim([0, 1]);
% xlabel('Time (Number of cycles)', 'FontName', 'Times New Roman', ...
% 'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
% ylabel('First-order index', 'FontName', 'Times New Roman', ...
% 'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
% legend({'Single loop Algorithm 1', 'Single loop Algorithm 2'},'location', 'northwest', 'FontName', 'Times New Roman', ...
% 'FontSize',10,'Interpreter', 'LaTeX')
% title('$K_c$','FontName', 'Times New Roman', ...
% 'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
% set(gcf, 'Position',  [200, 300, 900, 600])
% set(gcf, 'PaperPosition', [0 0 7 5]); %Position plot at left hand corner with width 5 and height 5.
% set(gcf, 'PaperSize', [7 5]); %Set the paper to have width 5 and height 5.
% saveas(gcf, fullfile(fname,'\Fig_UniMaxF_3_3'), 'pdf') %Save figure
% 
% 
% figure
% set(0,'DefaultLineLineWidth',2)
% % plot(First_order_double_loop, 'r+');
% % hold on;
% plot(t,SingleFirst_Index_Alg1(:,4), 'm-'); hold on;
% plot(t,SingleFirst_Index_Alg2(:,4), 'b--'); hold on;
% % plot(t,First_order_Saltelli(2:end,1), 'c.-.'); hold on;
% % plot(t,First_order_FAST(:,4), 'r-.'); hold on;
% % plot(t,First_order_Sobol(2:end,4), 'k:');
% ylim([0, 1]);
% xlabel('Time (Number of cycles)', 'FontName', 'Times New Roman', ...
% 'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
% ylabel('First-order index', 'FontName', 'Times New Roman', ...
% 'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
% legend({'Single loop Algorithm 1', 'Single loop Algorithm 2'},'location', 'northwest', 'FontName', 'Times New Roman', ...
% 'FontSize',10,'Interpreter', 'LaTeX')
% title('$a_0$','FontName', 'Times New Roman', ...
% 'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
% set(gcf, 'Position',  [200, 300, 900, 600])
% set(gcf, 'PaperPosition', [0 0 7 5]); %Position plot at left hand corner with width 5 and height 5.
% set(gcf, 'PaperSize', [7 5]); %Set the paper to have width 5 and height 5.
% saveas(gcf, fullfile(fname,'\Fig_UniMaxF_3_4'), 'pdf') %Save figure
% 
% 
% 
% 
% figure
% set(0,'DefaultLineLineWidth',2)
% % plot(First_order_double_loop, 'r+');
% % hold on;
% plot(t,SingleFirst_Index_Alg1(:,5), 'm-'); hold on;
% plot(t,SingleFirst_Index_Alg2(:,5), 'b--'); hold on;
% % plot(t,First_order_Saltelli(2:end,1), 'c.-.'); hold on;
% % plot(t,First_order_FAST(:,5), 'r-.'); hold on;
% % plot(t,First_order_Sobol(2:end,5), 'k:');
% ylim([0, 1]);
% xlabel('Time (Number of cycles)', 'FontName', 'Times New Roman', ...
% 'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
% ylabel('First-order index', 'FontName', 'Times New Roman', ...
% 'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
% legend({'Single loop Algorithm 1', 'Single loop Algorithm 2'},'location', 'northwest', 'FontName', 'Times New Roman', ...
% 'FontSize',10,'Interpreter', 'LaTeX')
% title('$\mu_{F_{max1}}$','FontName', 'Times New Roman', ...
% 'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
% set(gcf, 'Position',  [200, 300, 900, 600])
% set(gcf, 'PaperPosition', [0 0 7 5]); %Position plot at left hand corner with width 5 and height 5.
% set(gcf, 'PaperSize', [7 5]); %Set the paper to have width 5 and height 5.
% saveas(gcf, fullfile(fname,'\Fig_UniMaxF_3_5'), 'pdf') %Save figure
% 
% 
% figure
% set(0,'DefaultLineLineWidth',2)
% % plot(First_order_double_loop, 'r+');
% % hold on;
% plot(t,SingleFirst_Index_Alg1(:,6), 'm-'); hold on;
% plot(t,SingleFirst_Index_Alg2(:,6), 'b--'); hold on;
% % plot(t,First_order_Saltelli(2:end,1), 'c.-.'); hold on;
% % plot(t,First_order_FAST(:,6), 'r-.'); hold on;
% % plot(t,First_order_Sobol(2:end,5), 'k:');
% ylim([0, 1]);
% xlabel('Time (Number of cycles)', 'FontName', 'Times New Roman', ...
% 'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
% ylabel('First-order index', 'FontName', 'Times New Roman', ...
% 'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
% legend({'Single loop Algorithm 1', 'Single loop Algorithm 2'},'location', 'northwest', 'FontName', 'Times New Roman', ...
% 'FontSize',10,'Interpreter', 'LaTeX')
% title('$\mu_{F_{max2}}$','FontName', 'Times New Roman', ...
% 'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
% set(gcf, 'Position',  [200, 300, 900, 600])
% set(gcf, 'PaperPosition', [0 0 7 5]); %Position plot at left hand corner with width 5 and height 5.
% set(gcf, 'PaperSize', [7 5]); %Set the paper to have width 5 and height 5.
% saveas(gcf, fullfile(fname,'\Fig_UniMaxF_3_6'), 'pdf') %Save figure
% 
% 
% figure
% set(0,'DefaultLineLineWidth',2)
% % plot(First_order_double_loop, 'r+');
% % hold on;
% plot(t,SingleFirst_Index_Alg1(:,7), 'm-'); hold on;
% plot(t,SingleFirst_Index_Alg2(:,7), 'b--'); hold on;
% % plot(t,First_order_Saltelli(2:end,1), 'c.-.'); hold on;
% % plot(t,First_order_FAST(:,7), 'r-.'); hold on;
% % plot(t,First_order_Sobol(2:end,5), 'k:');
% ylim([0, 1]);
% xlabel('Time (Number of cycles)', 'FontName', 'Times New Roman', ...
% 'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
% ylabel('First-order index', 'FontName', 'Times New Roman', ...
% 'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
% legend({'Single loop Algorithm 1', 'Single loop Algorithm 2'},'location', 'northwest', 'FontName', 'Times New Roman', ...
% 'FontSize',10,'Interpreter', 'LaTeX')
% title('$\mu_{F_{max3}}$','FontName', 'Times New Roman', ...
% 'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
% set(gcf, 'Position',  [200, 300, 900, 600])
% set(gcf, 'PaperPosition', [0 0 7 5]); %Position plot at left hand corner with width 5 and height 5.
% set(gcf, 'PaperSize', [7 5]); %Set the paper to have width 5 and height 5.
% saveas(gcf, fullfile(fname,'\Fig_UniMaxF_3_7'), 'pdf') %Save figure
% 
% 
% figure
% set(0,'DefaultLineLineWidth',2)
% % plot(First_order_double_loop, 'r+');
% % hold on;
% plot(t,SingleFirst_Index_Alg1(:,8), 'm-'); hold on;
% plot(t,SingleFirst_Index_Alg2(:,8), 'b--'); hold on;
% % plot(t,First_order_Saltelli(2:end,1), 'c.-.'); hold on;
% % plot(t,First_order_FAST(:,8), 'r-.'); hold on;
% % plot(t,First_order_Sobol(2:end,7), 'k:');
% ylim([0, 1]);
% xlabel('Time (Number of cycles)', 'FontName', 'Times New Roman', ...
% 'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
% ylabel('First-order index', 'FontName', 'Times New Roman', ...
% 'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
% legend({'Single loop Algorithm 1', 'Single loop Algorithm 2'},'location', 'northwest', 'FontName', 'Times New Roman', ...
% 'FontSize',10,'Interpreter', 'LaTeX')
% title('$\sigma_{F_{max1}}$','FontName', 'Times New Roman', ...
% 'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
% set(gcf, 'Position',  [200, 300, 900, 600])
% set(gcf, 'PaperPosition', [0 0 7 5]); %Position plot at left hand corner with width 5 and height 5.
% set(gcf, 'PaperSize', [7 5]); %Set the paper to have width 5 and height 5.
% saveas(gcf, fullfile(fname,'\Fig_UniMaxF_3_8'), 'pdf') %Save figure
% 
% figure
% set(0,'DefaultLineLineWidth',2)
% % plot(First_order_double_loop, 'r+');
% % hold on;
% plot(t,SingleFirst_Index_Alg1(:,9), 'm-'); hold on;
% plot(t,SingleFirst_Index_Alg2(:,9), 'b--'); hold on;
% % plot(t,First_order_Saltelli(2:end,1), 'c.-.'); hold on;
% % plot(t,First_order_FAST(:,9), 'r-.'); hold on;
% % plot(t,First_order_Sobol(2:end,7), 'k:');
% ylim([0, 1]);
% xlabel('Time (Number of cycles)', 'FontName', 'Times New Roman', ...
% 'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
% ylabel('First-order index', 'FontName', 'Times New Roman', ...
% 'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
% legend({'Single loop Algorithm 1', 'Single loop Algorithm 2'},'location', 'northwest', 'FontName', 'Times New Roman', ...
% 'FontSize',10,'Interpreter', 'LaTeX')
% title('$\sigma_{F_{max2}}$','FontName', 'Times New Roman', ...
% 'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
% set(gcf, 'Position',  [200, 300, 900, 600])
% set(gcf, 'PaperPosition', [0 0 7 5]); %Position plot at left hand corner with width 5 and height 5.
% set(gcf, 'PaperSize', [7 5]); %Set the paper to have width 5 and height 5.
% saveas(gcf, fullfile(fname,'\Fig_UniMaxF_3_9'), 'pdf') %Save figure
% 
% figure
% set(0,'DefaultLineLineWidth',2)
% % plot(First_order_double_loop, 'r+');
% % hold on;
% plot(t,SingleFirst_Index_Alg1(:,10), 'm-'); hold on;
% plot(t,SingleFirst_Index_Alg2(:,10), 'b--'); hold on;
% % plot(t,First_order_Saltelli(2:end,1), 'c.-.'); hold on;
% % plot(t,First_order_FAST(:,10), 'r-.'); hold on;
% % plot(t,First_order_Sobol(2:end,7), 'k:');
% ylim([0, 1]);
% xlabel('Time (Number of cycles)', 'FontName', 'Times New Roman', ...
% 'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
% ylabel('First-order index', 'FontName', 'Times New Roman', ...
% 'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
% legend({'Single loop Algorithm 1', 'Single loop Algorithm 2'},'location', 'northwest', 'FontName', 'Times New Roman', ...
% 'FontSize',10,'Interpreter', 'LaTeX')
% title('$\sigma_{F_{max3}}$','FontName', 'Times New Roman', ...
% 'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
% set(gcf, 'Position',  [200, 300, 900, 600])
% set(gcf, 'PaperPosition', [0 0 7 5]); %Position plot at left hand corner with width 5 and height 5.
% set(gcf, 'PaperSize', [7 5]); %Set the paper to have width 5 and height 5.
% saveas(gcf, fullfile(fname,'\Fig_UniMaxF_3_10'), 'pdf') %Save figure
% 
% 
% 
% figure
% set(0,'DefaultLineLineWidth',2)
% % plot(First_order_double_loop, 'r+');
% % hold on;
% plot(t,SingleFirst_Index_Alg1(:,11), 'm-'); hold on;
% plot(t,SingleFirst_Index_Alg2(:,11), 'b--'); hold on;
% % plot(t,First_order_Saltelli(2:end,1), 'c.-.'); hold on;
% % plot(t,First_order_FAST(:,11), 'r-.'); hold on;
% % plot(t,First_order_Sobol(2:end,9), 'k:');
% ylim([0, 1]);
% xlabel('Time (Number of cycles)', 'FontName', 'Times New Roman', ...
% 'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
% ylabel('First-order index', 'FontName', 'Times New Roman', ...
% 'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
% legend({'Single loop Algorithm 1', 'Single loop Algorithm 2'},'location', 'northwest', 'FontName', 'Times New Roman', ...
% 'FontSize',10,'Interpreter', 'LaTeX')
% title('$U_{F_{max}}$','FontName', 'Times New Roman', ...
% 'FontSize',16,'Color','k', 'Interpreter', 'LaTeX')
% set(gcf, 'Position',  [200, 300, 900, 600])
% set(gcf, 'PaperPosition', [0 0 7 5]); %Position plot at left hand corner with width 5 and height 5.
% set(gcf, 'PaperSize', [7 5]); %Set the paper to have width 5 and height 5.
% saveas(gcf, fullfile(fname,'\Fig_UniMaxF_3_11'), 'pdf') %Save figure
% 
% 
