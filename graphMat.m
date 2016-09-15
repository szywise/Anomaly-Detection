function [IncMat,B,Bus,Nbr,SortDegIdx,iRefBus] = graphMat(casedata)
% graphMat:	calculate incidence matrix, weighted Laplacian matrix,
%			Bus-struct, Nbr from the graph given by casedata. randn_xMat is
%			used when calculating the mean of M. 
define_constants;
global DEG M_NOMINAL
nBus = size(casedata.bus,1);
nBranch = size(casedata.branch,1);

IncMat = zeros(nBus,nBranch); % incident matrix
Bus = zeros(nBus,3); % col 1: DEG; col 2: M_NOMINAL; col 3: IF_OBS
Nbr = zeros(nBus,nBus); % indices of neighboring buses
for iBranch = 1:nBranch % Node: No, deg, <r_j>, M
	f_bus = casedata.branch(iBranch,F_BUS);
	t_bus = casedata.branch(iBranch,T_BUS);
	IncMat(f_bus,iBranch) = 1;
	IncMat(t_bus,iBranch) = -1;	
	Bus(f_bus,DEG) = Bus(f_bus,DEG) + 1;
	Nbr(f_bus,Bus(f_bus,DEG)) = t_bus;
	Bus(t_bus,DEG) = Bus(t_bus,DEG) + 1;
	Nbr(t_bus,Bus(t_bus,DEG)) = f_bus;
end

Bus(:,M_NOMINAL) = calcMetric(casedata,zeros(nBranch,1),0); % nominal of M

[~,SortDegIdx] = sort(Bus(:,DEG),'descend'); % if [B I]=sort(A), then A[I]=B.
iRefBus = SortDegIdx(1); % index of the reference bus
IncMat(iRefBus,:) = []; % make 1st observable bus the reference bus, and remove the row

DiagMat = diag(1 ./ casedata.branch(:,BR_X)); % nBr-by-nBr diagnal matrix
B = IncMat*DiagMat*IncMat'; % Weighted Laplacian matrix, (N-1)by(N-1)