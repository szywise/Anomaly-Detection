function pos = reconstr(B,IncMat,Obs,diff_theta)
% Sparse coefficient reconstruction using SVD and OMP

invB = inv(B); % (N-1)by(N-1)
iRefBus = Obs(1);
Obs_B = Obs(2:end); % observed bus indices without reference bus
Obs_B(Obs_B>iRefBus) = Obs_B(Obs_B>iRefBus)-1; % adjust indices
[U,S,V] = svd(invB(Obs_B,:),'econ');
y = S\U'*diff_theta(Obs(2:end)); % measurement vector
A = V'*IncMat; % dictionary, (I-1)by(L)

normA = sqrt(sum(A.*A))'; % norm of each column of A
res = y; % residual error

[~,max_idx] = max(abs(A'*res)./normA); % find max's index
pos = max_idx; % store the index

