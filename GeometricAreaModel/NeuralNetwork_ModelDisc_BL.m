

sweep = [210;28;1;1];                 % parameter values to test
scores = zeros(length(sweep), 1);       % pre-allocation
models = cell(length(sweep), 1);        % pre-allocation
x = Inp(:,1:64);                             % inputs
t = Out(:,1:64);                             % targets
trainFcn = 'trainscg';                  % scaled conjugate gradient
for i = 1:length(sweep)
    hiddenLayerSize = sweep(i);         % number of hidden layer neurons
    net = patternnet(hiddenLayerSize);  % pattern recognition network
    net.divideParam.trainRatio = 70/100;% 70% of data for training
    net.divideParam.valRatio = 15/100;  % 15% of data for validation
    net.divideParam.testRatio = 15/100; % 15% of data for testing
    net = train(net, x, t);             % train the network
    models{i} = net;                    % store the trained network
    p = net(Inp(:,65:end));                     % predictions
    [~, p] = max(p);                    % predicted labels
    scores(i) = sum(Out(:,65:end) == p) /...    % categorization accuracy
        length(Out(:,65:end));
end