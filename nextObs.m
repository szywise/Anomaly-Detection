function [Buff] = nextObs(nPer,diff_theta,Nbr,Repo)
% nextObs:		choose buses for next measurement
% input param:	nPer, diff_theta, Bus, Neighbors, (no Obs), Repo
% output param:	Buff
global DEG M_NOMINAL IF_OBS
global Bus
Buff = [];
while length(Buff)<nPer
	[~, jj] = max(abs(diff_theta(Repo))); % max deviation
	ii = Repo(jj); % bus index with largest phase angle deviation
	NbrList = Nbr(ii,1:Bus(ii,DEG)); % Neighbor list
	NbrList = NbrList(Bus(NbrList,IF_OBS)==0); % unobserved neighbors, better than setdiff
	NbrList = unique(NbrList); % there are two lines connecting 89 and 92
	
	[~,SortMetricIdx]=sort(Bus(NbrList,M_NOMINAL),'descend'); % sort based on decreasing M
	dBuff = NbrList(SortMetricIdx); % neighboring buses' indices ranked by Metric
	if isempty(dBuff)
		dBuff = zeros(1,0); % eliminate warning of row inconsistency of concatenated array
	end
	if length(Buff) + length(dBuff) < nPer
		Buff = [Buff dBuff];
		Repo(jj) = []; % remove expanded(used) bus from Repo
	else
		Buff = [Buff dBuff(1:nPer-length(Buff))]; % not fully expanded, keep it
	end
	Bus(Buff,IF_OBS) = 1; % tag as 'observed'
end