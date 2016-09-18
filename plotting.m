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
perc = [1 5 10 15];
for ii = 1:4
	subplot(2,2,ii); stem(diff_theta_order(:,ii));
	title([num2str(perc(ii)) '% reactance fluctuation']);
	ylim([-15 15]); grid on; 
	xlabel('Bus No.'); ylabel({'Changes of |\Delta\theta| ranking','compared to 1% fluctuation'});
end
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