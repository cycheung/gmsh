     
A = importdata('temp.txt');

subplot(2,2,[1 2])
plot( A(:,1), A(:,2), 'k', 'Linewidth', 1.5); hold on
plot( A(:,1), A(:,8), 'k', 'Linewidth', 1.5); hold on
plot( A(:,1), A(:,14), 'k', 'Linewidth', 1.5); hold on
plot( A(:,1), A(:,20), 'k', 'Linewidth', 1.5); hold on
plot( A(:,1), A(:,26), 'k', 'Linewidth', 1.5); hold on
plot( A(:,1), 320.15*ones(size(A(:,1))), 'LineStyle', '-', 'Color', [0.8 0.8 0.8], 'Linewidth', 1); hold on

legend('z = 0 mm', 'z = 0.05 mm', 'z = 0.10 mm', 'z = 0.15 mm', 'z = 0.20 mm', 'Seuil fibres A\delta', 'Seuil fibres C');
xlabel('Temps  [s]'); ylabel('Température  [K]');
title('Température de la peau exposée à un laser CO2 pour différentes profondeurs');


