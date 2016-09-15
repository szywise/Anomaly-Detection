function [Metric] = calcMetric(casedata, noisevec, xSigMax)
% Calculate metric from casedata with noise
BR_X = 4;
F_BUS = 1;
T_BUS = 2;
nBus = size(casedata.bus,1);
nBranch = size(casedata.branch,1);
xRange = (2*xSigMax).*noisevec + repmat(1-xSigMax,nBranch,1); % (1-xSig,1+xSig)
awgnBrX = casedata.branch(:,BR_X) .* xRange; % branch reactance with Gaussian noise
Degree = zeros(nBus,1);
Metric = zeros(nBus,1);
Nbr = zeros(nBus,nBus); % reactance of neighbors
for jBranch = 1:nBranch
	f_bus = casedata.branch(jBranch,F_BUS);
	t_bus = casedata.branch(jBranch,T_BUS);
	Degree(f_bus) = Degree(f_bus) + 1;
	Nbr(f_bus,Degree(f_bus)) = awgnBrX(jBranch);
	Degree(t_bus) = Degree(t_bus) + 1;
	Nbr(t_bus,Degree(t_bus)) = awgnBrX(jBranch);
end
for iBus = 1:nBus
	if Degree(iBus) > 1
		XNbr = Nbr(iBus,1:Degree(iBus)); % reactance of neighboring lines
		beta = 1/sum(1./XNbr); % beta-value of the bus
		RNbr = beta./XNbr; % r-value of neighboring lines
		Metric(iBus) = mean(-log(1-RNbr.^2));
	end
end