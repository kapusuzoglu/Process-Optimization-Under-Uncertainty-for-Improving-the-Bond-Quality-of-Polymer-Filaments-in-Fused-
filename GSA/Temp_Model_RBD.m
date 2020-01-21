function [T2] = ...
         Temp_Model_RBD(input, other_inputs)
     
Mu = other_inputs(1:4);
Std = other_inputs(5:8);

Other_inputs=other_inputs(1,9:15);

rho = norminv(input(:, 1), Mu(1), Std(1));
C = norminv(input(:, 2), Mu(2), Std(2));
hc = norminv(input(:, 3), Mu(3), Std(3));
hk = norminv(input(:, 4), Mu(4), Std(4));

Mat = [rho,C,hc,hk];

N_Samples=size(input,1);
parfor ii=1:N_Samples
    T2(ii,:) = FDM_TempModel_GSA_fil2(Other_inputs, Mat(ii,:) );
end

end