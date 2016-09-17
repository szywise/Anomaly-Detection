function pos = reconstr(B,B_true,IncMat,Obs,diff_theta)
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

invB1 = inv(B_true);
[U1,S1,V1] = svd(invB1(Obs_B,:),'econ');
y1 = S1\U1'*diff_theta(Obs(2:end));
% A1 = A*6;
A1 = V1'*IncMat;
normA1 = sqrt(sum(A1.*A1))';

A_corr = acos(diag(A'*A1) ./ (normA1.*normA))*180/pi;
y_corr = acos(y1'*y/(norm(y1)*norm(y)))*180/pi

[~,max_idx] = max(abs(A'*res)./normA); % find max's index
pos = max_idx; % store the index

