function index = SingleLoop_FirstOrder_Alg1( input, output, ndomain, varargin)
% input: nsample * nd matrix, where nsample is the number of sample, and nd
% is the input dimension
% output: nsample * 1 array
% ndomain: number of sub-domain to divide a single input, sqrt(nsample)
% recommended and nsample larger than 50 recommended.
% varargin: incudes two arg. arg 1 is a vector of dimension number of
% inputs for which cdf functions are available. For example, [2, 3] means
% input 2 and input 3 have cdf functions. arg 2 is a cell of corresponding cdf functions

% This algorithm is proposed by me. Please cite the following paper if you use this code.
% Li, Chenzhao, and Sankaran Mahadevan. "An efficient modularized sample-based method to estimate the first-order Sobol’index." Reliability Engineering & System Safety (2016).


[nsample, nd] = size(input);

%% convert the samples into cdf domains
U = linspace(0, 1, ndomain+1);
U = U';

cdf_input = zeros(nsample, nd);

cdf_values = transpose(linspace(1/nsample, 1, nsample));
if nargin == 3 % cdf function provided
    varargin = {[], {}};
end

j = 1;
ismember(1, varargin{1}{1});

for i = 1:nd
    if ismember(i, varargin{1}{1})
        cdf_input(:, i) = cdf(varargin{1}{2}{j}, input(:, i), varargin{1}{2+i}{1}, varargin{1}{2+i}{2});
        j = j+1;
    else
        [~, IX] = sort(input(:, i));
        [~, IX2] = sort(IX);
        cdf_input(:, i) = cdf_values(IX2);
    end
end

% unique(cdf_input(:,2))
% figure(3)
% plot(unique(cdf_input(:,1)))


%% compute indices
VY = var(output);
MeanY_local = zeros(ndomain, nd, size(output,2));
for i = 1:nd
    cdf_input_i = cdf_input(:, i);
    output_i = output;
    U_i = U;
    for j = 1:ndomain
        sub = cdf_input_i<U_i(j+1);
        MeanY_local(j, i, :) = mean(output_i(sub,:));
        inverse_sub = ~sub;
        cdf_input_i = cdf_input_i(inverse_sub);
        output_i = output_i(inverse_sub,:);
    end
end
% index = mean(VarY_local)/VY;

%index = var(MeanY_local)./VY;
index = zeros(size(output,2)-1, nd);
MeanY_local = MeanY_local(:,:,2:end);
VY = VY(2:end);

% % meanY = MeanY_local(:,:,1)
% % meanY5 = MeanY_local(:,:,5000)

% var(MeanY_local(:,:,1))
% VY(1,1)
% % VY(1,5000)
for jj = 1:size(output,2)-1
    index(jj,:) = var(MeanY_local(:,:,jj))/VY(1,jj);
end


end