function [s,pos,res] = reconstr(B,IncMat,Obs,diff_theta,kappa)
% Sparse coefficient reconstruction using SVD and OMP
nBranch = size(IncMat,2); % set it global?
invB = inv(B); % (N-1)by(N-1)
iRefBus = Obs(1);
Obs_B = Obs(2:end); % observed bus indices without reference bus
Obs_B(Obs_B>iRefBus) = Obs_B(Obs_B>iRefBus)-1; % adjust indices
[U,S,V] = svd(invB(Obs_B,:),'econ');
y = S\U'*diff_theta(Obs(2:end)); % measurement vector
A = V'*IncMat; % dictionary, (I-1)by(L)

normA = sqrt(sum(A.*A))'; % norm of each column of A
res = y; % residual error
pos = zeros(1,kappa); % position: the indices of nonzeros entries of s
while 1 % norm(res)>gamma ... I don't understand
	for kk = 1:kappa % for kk-th iteration
		[~,max_idx] = max(abs(A'*res)./normA); % find max's index
		pos(kk) = max_idx; % store the index
		Ak = A(:,pos(1:kk)); % chosen atoms from dictionary
		s = zeros(nBranch,1); % initialize coefficient s
		s(pos(1:kk)) = (Ak'*Ak)\Ak'*y; % least-square
		res = y - A*s; % residual error
	end
	break;
end