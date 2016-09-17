close
xfluc = [0 5 10 15 20 25 30 35 40 45 60 80 90];
acc = [0.7944 0.79389 0.78922 0.79215 0.79349 0.79320 0.78982 0.79054 0.79252 0.79144 0.788 0.76242 0.74498];
plot(xfluc, acc, 'o--');
ylim([0.75 0.80]);
grid on;
xlabel('reactance fluctuation (%)');
ylabel('accuracy');
