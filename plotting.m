%% accuracy respond to reactance fluctuation level
% close
% xfluc = [0 5 10 15 20 25 30 35 40 45 60 80 90];
% acc = [0.7944 0.79389 0.78922 0.79215 0.79349 0.79320 0.78982 0.79054 0.79252 0.79144 0.788 0.76242 0.74498];
% plot(xfluc, acc, 'o--');
% ylim([0.75 0.80]);
% grid on;
% xlabel('reactance fluctuation (%)');
% ylabel('accuracy');
%% change of M ranking with different reactance fluctuation level
% perc = [1 5 10 15];
% for ii = 1:4
% 	subplot(2,2,ii); stem(diff_ranking(ii,:));
% 	title([num2str(perc(ii)) '% reactance fluctuation']);
% 	ylim([-15 15]); grid on; 
% 	xlabel('Bus No.'); ylabel('Changes of M ranking');
% 	text(60,12,['average ranking deviation: ',num2str(mean(abs(diff_ranking(ii,:))))],...
% 		'HorizontalAlignment','center');
% end
%% change of diff_theta ranking with different reactance fluctuation level
% perc = [1 5 10 15];
% for ii = 1:4
% 	subplot(2,2,ii); stem(diff_theta_order(:,ii));
% 	title([num2str(perc(ii)) '% reactance fluctuation']);
% 	ylim([-15 15]); grid on; 
% 	xlabel('Bus No.'); ylabel({'Changes of |\Delta\theta| ranking','compared to 1% fluctuation'});
% end
%% diff_theta of different reactance fluctuation level
% for ii = [10 5 8 9]%1:n_xSig
% 	plot(rec_diff_theta(ii,:),'+-');
% 	grid on; hold on;
% end
% title('Phasor angle deviation');
% xlabel('Bus No.');ylabel('\Delta\theta','interpreter','tex');
% legend('0','0.5','0.8','0.9');%,'0.5','0.6','0.7','0.8','0.9'
%% change of M ranking with different reactance fluctuation level
% perc = [1 5 10 15];
% for ii = 1:4
% % 	subplot(2,2,ii); 
% 	plot(rec_diff_theta(ii,:),'+-');
% 	title([num2str(perc(ii)) '% reactance fluctuation']);
% 	ylim([-0.12 0.12]); grid on; 
% 	xlabel('Bus No.'); ylabel('phasor angle deviation');
% 	hold on;
% end
%% sorted nominal M and corresponding actual M
% plot(M(M_order),'+'); hold on; plot(M1(M_order),'d-');legend('sorted nominal M','corresponding actual M')
% title('M changes with 15% reactance fluctuation')
% grid on;
%% Laplacian matrix
% figure(1);
% subplot(221); imagesc(pre_Lap); colorbar; title('pre_Lap');
% subplot(222); imagesc(post_Lap); colorbar; title('post_Lap');
% subplot(223); imagesc(inv(pre_Lap)); colorbar; title('inv pre_Lap');
% subplot(224); imagesc(inv(post_Lap)); colorbar; title('inv post_Lap');
% figure(2);
% subplot(121);
% imagesc(inv(pre_Lap)); colorbar;
% title({'(a) Inverse of Laplacian matrix';'B^{-1}'});
% subplot(122);
% dinvLap = inv(pre_Lap)-inv(post_Lap);
% imagesc(dinvLap,[-.03 .02]); colorbar;
% title({'(b) Difference of inverse of Laplacian matrix';'B^{\prime -1}-B^{-1}'});
% figure(3);
% pre_power(iRefBus) = [];
% imagesc(dinvLap*pre_power); colorbar; title('\Delta\theta');
% figure(4);
% imagesc(inv(pre_Lap),[0,0.4]); colormap parula; colorbar;
% title('Inverse of Laplacian Matrix');
% figure(5);
% inv_pre_Lap = inv(pre_Lap);
% plot(sort(inv_pre_Lap(:)),'*-');
%% distribution of Laplacian matrix
% plot(sort(inv_pre_Lap(:)),'*-');
% grid on;
% title('distribution of inv(B)');
%% Model verification
close all;
plot(diff_theta); hold on;

H = inv(pre_Lap);
m = IncMat(:,31);
p = mpc_init.bus(:,PD);
NoRef = [1:iRefBus-1, iRefBus+1:118];
c = mpc_init.branch(31,BR_X)-m'*H*m;
diff_theta_theory = zeros(118,1);
diff_theta_theory(NoRef) = H*m*m'*H*p(NoRef)./c;
diff_theta_theory(iRefBus) = diff_theta(iRefBus);
plot(diff_theta_theory); hold on;

plot(diff_theta_theory-3244/1559*diff_theta);
legend('experimental','th1','th2');
grid on; title('phasor angle deviation');
%% bug narrow down
plot(d