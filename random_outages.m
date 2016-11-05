% FUNCTION:	Performance demanding
% USAGE:	nAll	: number of measurements
%			xSigMax	: reactance fluctuation level
% MODIFIED:	Sep 18
clear
define_constants;
global DEG M_NOMINAL IF_OBS
global Bus
DEG = 1; % degree
M_NOMINAL = 2; % nominal value of M
IF_OBS = 3; % if-observed tag

%% Operation perameters
% make all parameters random
nAll = 50; % number of meansurements
nPer = 5; % number of measurements each time instant
arr_xSig = 0; %[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0]; % reactance fluctuation level
n_xSig = length(arr_xSig);
n_xNoi = 1; %10; % realization of noise
n_pNoi = 1;
pSig = 0.0; % power noise level
load rand_x;
load randn_p;

%% Start
mpc_init = loadcase('case118');
nBus = size(mpc_init.bus,1); % number of buses(node)
nBranch = size(mpc_init.branch,1); % number of branches(edge)
Unsucc_Br = [7,9,113,133,134,176,177,183,184];
Succ_Br = setdiff(1:nBranch,Unsucc_Br);

%% Graph matrices: incidence matrix, Bus(DEG,M_MEAN), and Nbr
[IncMat,B,Bus,Nbr,Buff_init,iRefBus] = graphMat(mpc_init,nPer);
[~,~,M_ranking] = unique(Bus(:,M_NOMINAL));

%% Single line outage
arr_pAcc = zeros(size(arr_xSig));

rec_ranking_M = zeros(n_xSig,nBus);
rec_change_ranking_M = zeros(n_xSig,nBus);

rec_diff_theta = zeros(n_xSig,nBus);
rec_ranking_diff_theta = zeros(n_xSig,nBus);

rec_Obs = zeros(n_xSig,nAll);
rec_pos = zeros(2,nAll/nPer,n_xSig);

for i_xSig = 1:length(arr_xSig)
	start_time = datestr(now,'mmm dd, HH:MM:SS');
	cntAll = 0; % accumulate during iteration
	cntAcc = 0; % count of good recovery
	xSig = arr_xSig(i_xSig);
	for i_xNoi = 1:n_xNoi
		mpc_x = mpc_init; % backtracking needs more time
rand_x = rand(nBranch,1);
% rand_x = randi([0,1],nBranch,1);
		xRange = (2*xSig).*rand_x + repmat((1-xSig),nBranch,1); % (1-xSig,1+xSig)
		mpc_x.branch(:,BR_X) = mpc_x.branch(:,BR_X) .* xRange;
pre_Lap = IncMat*diag(1./mpc_x.branch(:,BR_X))*IncMat';

M1 = calcMetric(mpc_x,zeros(nBranch,1),0);
B1 = IncMat*diag(mpc_x.branch(:,BR_X))*IncMat';
[~,~,rec_ranking_M(i_xSig,:)] = unique(M1');
rec_change_ranking_M(i_xSig,:) = rec_ranking_M(i_xSig,:) - M_ranking';
[xSig norm(rec_change_ranking_M)/sqrt(118) mean(abs(rec_change_ranking_M))];

		[mpr_x,~] = rundcpf(mpc_x,mpoption('verbose',0,'out.all',0));
		pre_theta = mpr_x.bus(:,VA) * pi / 180; % all phasor angles
% 		pre_theta = pre_theta - repmat(pre_theta(iRefBus),nBus,1);
		pre_power = mpr_x.bus(:,PD);
		
		arr_Branch = 31;%Succ_Br(randperm(length(Succ_Br),20)); % randperm > datasample; no sort.
		for iBranch = arr_Branch
			mpc_xp = mpc_x;
diff_Lap = IncMat(:,iBranch) * (1./mpc_xp.branch(iBranch,BR_X)) * IncMat(:,iBranch)';
post_Lap = pre_Lap - diff_Lap;
			mpc_xp.branch(iBranch,:) = []; % branch iBr breaks down
			for i_pNoi = 1:n_pNoi
randn_p = randn(nBus,1);
				pNoise = pSig * mean(pre_power) * randn_p; % randn(nBus,1);
				mpc_xp.bus(:,PD) = pre_power + pNoise;

				% run optimal power flow simulation
				[mpr_xp,succ] = rundcpf(mpc_xp,mpoption('verbose',0,'out.all',0));
				if ~succ
					break; % this topology is disconnected, go to next
				end
				post_theta = mpr_xp.bus(:,VA) * pi / 180;
% 				post_theta = post_theta - repmat(post_theta(iRefBus),nBus,1);
				diff_theta = post_theta - pre_theta;
rec_diff_theta(i_xSig,:) = diff_theta';
[~,~,rec_ranking_diff_theta(i_xSig,:)] = unique(abs(diff_theta'));
rec_ranking_diff_theta(i_xSig,:) = nBus+1-rec_ranking_diff_theta(i_xSig,:);
				% data-acquisition and decision-making
				Bus(:,IF_OBS) = 0;
				Buff = Buff_init; % init dynamic measurement buffer
				Bus(Buff,IF_OBS) = 1; % tag these buses as 'observed'
				Obs = []; % observed buses
				Repo = []; % observed but not expanded buses, just like a repository    
				cntTime = 0;
				while 1
					cntTime = cntTime + 1;
					Obs = [Obs Buff];
					Repo = [Repo Buff];

					% sparse coefficient reconstruction: SVD & OMP
					pos = reconstr(B,IncMat,Obs,diff_theta); 
pos1 = reconstr(B1, IncMat,Obs,diff_theta);
rec_pos(:,cntTime,i_xSig) = [pos;pos1];

					if length(Obs) >= nAll % number of measurements are fixed
						break; % 100 internal buses
					end

					% choose buses for next measurement
					Buff = nextObs(nPer,diff_theta,Nbr,Repo);
				end % end of data-acquisition and decision-making

				cntAll = cntAll + 1;
				if pos == iBranch
					cntAcc = cntAcc +  1;
				end
			end
% 			[iBranch pos(1)]
			rec_Obs(i_xSig,:) = Obs;
		end
		[i_xSig i_xNoi];
	end % for each branch
	pAcc = cntAcc / cntAll 
	arr_pAcc(i_xSig) = pAcc;
	%% Output file header
	fid = fopen('accuracy_log.txt','a');
	fprintf(fid, ['Start: ',start_time,' - Stop: ',datestr(now,'mmm dd, HH:MM:SS'),'\n']);
	fprintf(fid, 'Measurements: %d\tn_xNoi: %d\tn_pNoi: %d\n',nAll,n_xNoi,n_pNoi);
	fprintf(fid, 'xSig: %d%%\t\t\tAccuracy: %f%%\n',xSig*100,pAcc*100);
	fprintf(fid,'\n');
	fclose(fid);
end
rec_sort_Obs = sort(rec_Obs,2);
[~,sort_dt_idx]=sort(abs(diff_theta),'descend');
largest = sort_dt_idx(1:50)';
rec_concern = rec_ranking_diff_theta(:,largest);
% rec_concern_change_ranking = rec_concern - repmat(rec_concern(10,:),10,1);