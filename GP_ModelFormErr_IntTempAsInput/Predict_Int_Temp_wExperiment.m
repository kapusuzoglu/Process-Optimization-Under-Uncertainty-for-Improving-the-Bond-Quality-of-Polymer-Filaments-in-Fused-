clc;clear;
close all;

% ======================================================================= %
%__________________________       INPUTS     ______________________________
%Process Variables
layer = 6; % which layer we are looking at for the bond lengths 
% numFilaments-1: (numFilaments-1)th interface on that given layer
numFilaments = 15;
%__________________________________________________________________________

Sample_Temp_Speed = 1.0*[210 28; 210 28;...
                         234 34; 234 34;...
                         219 23;...
                         242 26;...
                         260 35; 260 35;...
                         223 29; 223 29;...
                         247 37; 247 37;...
                         255 46; 255 46; ...
                         252 21; 240 42; 240 42; 240 42;
                         239 29; 239 40.4; 239 36.6; 239 32.8; 239 25.2;
                         226 18.6; 226 18.6; 240 15];


NumTest = size(Sample_Temp_Speed,1);


KK = 273.15; % 0 Celcius in K
% convert Celcius to Kelvin

% ======================================================================= %
%% Experiment
%% This will plot the temperatures at line B (x = 175 pixels)
% Row pixels that are manually chosen which represents the interface 
% between lines and column pixels represent point B
row=[41,38,35,32,30,27,25,22,20,18,15,13,10,8];
freq = 10.5; % frequency of frames

% [t210s28_1_temps.mat; t210s28_2_temps.mat; ;
% t234s34_4_temps.mat; t234s34_2_temps.mat
% t219s23_5_temps.mat; 
% t242s26_3_temps.mat; 
% t260s35_1_temps.mat; t260s35_3_temps.mat;
% t223s29_5_temps.mat; t223s29_6_temps.mat
% t247s37_2_temps.mat; t247s37_3_temps.mat;
% t255s46_2_temps.mat; t255s46_3_temps.mat;
% t252s21_3_temps.mat;
% t240s42_1_temps.mat; t240s42_2_temps.mat; t240s42_3_temps.mat;
% t239s29_1_temps.mat
% t239s40.4_1; t239s36.6_1; t239s32.8_1; t239s25.2_1;
% t226s18_6_1; t226s18_6_2; t240s15_1]
Temp.S1 = load('..\Experiment_Analytical_Abaqus\T210S28-selected\t210s28_1_temps.mat');
Temp.S2 = load('..\Experiment_Analytical_Abaqus\T210S28-selected\t210s28_2_temps.mat');
Temp.S3 = load('..\Experiment_Analytical_Abaqus\T234S34-selected\t234s34_4_temps.mat');
Temp.S4 = load('..\Experiment_Analytical_Abaqus\T234S34-selected\t234s34_2_temps.mat');
Temp.S5 = load('..\Experiment_Analytical_Abaqus\T219S23-selected\t219s23_5_temps.mat');
Temp.S6 = load('..\Experiment_Analytical_Abaqus\T242S26-selected\t242s26_3_temps.mat');
Temp.S7 = load('..\Experiment_Analytical_Abaqus\T260S35-selected\t260s35_1_temps.mat');
Temp.S8 = load('..\Experiment_Analytical_Abaqus\T260S35-selected\t260s35_3_temps.mat');
Temp.S9 = load('..\Experiment_Analytical_Abaqus\T223S29-selected\t223s29_5_temps.mat');
Temp.S10 = load('..\Experiment_Analytical_Abaqus\T223S29-selected\t223s29_6_temps.mat');
Temp.S11 = load('..\Experiment_Analytical_Abaqus\T247S37-selected\t247s37_2_temps.mat');
Temp.S12 = load('..\Experiment_Analytical_Abaqus\T247S37-selected\t247s37_3_temps.mat');
Temp.S13 = load('..\Experiment_Analytical_Abaqus\T255S46-selected\t255s46_2_temps.mat');
Temp.S14 = load('..\Experiment_Analytical_Abaqus\T255S46-selected\t255s46_3_temps.mat');
Temp.S15 = load('..\Experiment_Analytical_Abaqus\T252S21-selected\t252s21_3_temps.mat');
Temp.S16 = load('..\Experiment_Analytical_Abaqus\T240S42\t240s42_1_temps.mat');
Temp.S17 = load('..\Experiment_Analytical_Abaqus\T240S42\t240s42_2_temps.mat');
Temp.S18 = load('..\Experiment_Analytical_Abaqus\T240S42\t240s42_3_temps.mat');
Temp.S19 = load('..\Experiment_Analytical_Abaqus\T239-selected\S29\t239s29_1_temps.mat');
Temp.S20 = load('..\Experiment_Analytical_Abaqus\T239-selected\S40.4\t239s40_4_1_temps.mat');
Temp.S21 = load('..\Experiment_Analytical_Abaqus\T239-selected\S36.6\t239s36_6_1_temps.mat');
Temp.S22 = load('..\Experiment_Analytical_Abaqus\T239-selected\S32.8\t239s32_8_1_temps.mat');
Temp.S23 = load('..\Experiment_Analytical_Abaqus\T239-selected\S25.2\t239s25_2_1_temps.mat');
Temp.S24 = load('..\Experiment_Analytical_Abaqus\T226S18.6\t226s18_6_1_temps.mat');
Temp.S25 = load('..\Experiment_Analytical_Abaqus\T226S18.6\t226s18_6_2_temps.mat');
Temp.S26 = load('..\Experiment_Analytical_Abaqus\S15-selected\T240\t240s15_1_temps.mat');


% inputs for GP model
t_final=zeros(layer*(numFilaments-1),NumTest);
slope=zeros(layer*(numFilaments-1),NumTest);
intercept=zeros(layer*(numFilaments-1),NumTest);


tfinal=zeros(numFilaments-1,NumTest);
slop=zeros(numFilaments-1,NumTest);
interc=zeros(numFilaments-1,NumTest);

First_idx=zeros(numFilaments-1,1);
Last_idx=zeros(numFilaments-1,1);

for iii=1:NumTest
    % import temperature data
    temps = Temp.(strcat('S', int2str(iii)));
    tic
    for jj=1:layer

        % Note: numFilaments-1 = size(row,2)
        pixeltemp = zeros(size(row,2), size(temps.(strcat('layer', int2str(jj))),3));
        Int_temps_exp = zeros(size(row,2), size(temps.(strcat('layer', int2str(jj))),3));
        FrameIdx=zeros(1,numFilaments-1);
        kk=1;
        for rr = 1:size(row,2)
            % convert Celcius to Kelvin, 174:180 -> find the max temp. near mid
            % point, col=175 is not always the mid point as we go to next lines bc. of
            % angle of the camera
            pixeltemp(rr, :) = max(temps.(strcat('layer', int2str(jj)))(row(rr), 174:180, :))+KK;
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
%             kkk=(jj-1)*(numFilaments-1)+rr;
            First_idx(rr) = find(Int_temps_exp(rr,:)~=0, 1, 'first');
            Last_idx(rr) = find(Int_temps_exp(rr,:)>130+KK, 1, 'last');
            
%             % remove zero elements
%             Int_temps_exp( :, ~any(Int_temps_exp,1) ) = [];  %columns
            
            
            start_frame = 1;
            final_frame = Last_idx(rr);
            time = start_frame/freq:1/freq:final_frame/freq;
            
            xx = time(First_idx(rr):Last_idx(rr));
            yy = Int_temps_exp(rr,First_idx(rr):Last_idx(rr));
%             xx=time;yy=Int_temps_exp(rr,:);
            [slop(rr,iii), interc(rr,iii)] = logfit(xx,yy,'logy');
            yApprox = exp(interc(rr,iii))*exp(slop(rr,iii)*xx);
            
            % save final time (which the filaments reach 130)
            tfinal(rr,iii) = xx(end);
        
        end
        kkk=(jj-1)*(numFilaments-1)+1; kkkk=(jj-1)*(numFilaments-1)+rr;
        t_final(kkk:kkkk,iii) = tfinal(:,iii);
        slope(kkk:kkkk,iii) = slop(:,iii);
        intercept(kkk:kkkk,iii) = interc(:,iii);
 
        
    end
    

    
toc
end


save('Experiment.mat','t_final','slope','intercept')



