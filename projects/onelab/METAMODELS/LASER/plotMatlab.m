
A = importdata('Tmax_00.txt');
B = importdata('Tmax_05.txt');
C = importdata('Tmax_10.txt');
D = importdata('Tmax_15.txt');
E = importdata('Tmax_20.txt');

dt = A(1,5);

subplot(2,2,[1 2])
plot([0; A(:,2)*dt+dt], [305.15; A(:,8)], 'k', 'Linewidth', 1.5); hold on
plot([0; B(:,2)*dt+dt], [305.15; B(:,8)], '--k', 'Linewidth', 1); hold on
plot([0; C(:,2)*dt+dt], [305.15; C(:,8)], ':k', 'Linewidth', 1); hold on
plot([0; D(:,2)*dt+dt], [305.15; D(:,8)], '-k', 'Linewidth', 1); hold on
plot([0; E(:,2)*dt+dt], [305.15; E(:,8)], '-.k', 'Linewidth', 1); hold on
plot(E(:,2)*dt, 320.15*ones(size(E(:,2))), 'LineStyle', '-', 'Color', [0.8 0.8 0.8], 'Linewidth', 1); hold on
plot(E(:,2)*dt, 314.15*ones(size(E(:,2))), 'LineStyle', '--', 'Color', [0.8 0.8 0.8], 'Linewidth', 1); hold on

legend('z = 0 mm', 'z = 0.05 mm', 'z = 0.10 mm', 'z = 0.15 mm', 'z = 0.20 mm', 'Seuil fibres A\delta', 'Seuil fibres C');
xlabel('Temps  [s]'); ylabel('Température  [K]');
title('Température de la peau exposée à un laser CO_2 pour différentes profondeurs');

F = importdata('Area_00_ad.txt');
G = importdata('Area_05_ad.txt');
H = importdata('Area_10_ad.txt');
I = importdata('Area_15_ad.txt');
J = importdata('Area_20_ad.txt');

Mat = [F(:,8)  G(:,8)  H(:,8)  I(:,8)  J(:,8)];
n = length(F);
l = [n n n n n];

for i=1:5
    for j=1:length(F)
        if (Mat(j,i)==0)
            l(i)=l(i)-1;
        end
    end
end

x = [0 0.05 0.10 0.15 0.20];
y1 = [max(F(:,8)) max(G(:,8)) max(H(:,8)) max(I(:,8)) max(J(:,8))];

K = importdata('Area_00_c.txt');
L = importdata('Area_05_c.txt');
M = importdata('Area_10_c.txt');
N = importdata('Area_15_c.txt');
O = importdata('Area_20_c.txt');

Mat = [K(:,8)  L(:,8)  M(:,8)  N(:,8)  O(:,8)];
n = length(K);
m = [n n n n n];

for i=1:5
    for j=1:length(K)
        if (Mat(j,i)==0)
            m(i)=m(i)-1;
        end
    end
end

y2 = [max(K(:,8)) max(L(:,8)) max(M(:,8)) max(N(:,8)) max(O(:,8))];

subplot(2,2,3)
plot(x, y1*10^6, 'k', x, y1*10^6, '.k', 'Linewidth', 1, 'MarkerSize', 14); hold on
plot(x, y2*10^6, 'k', x, y2*10^6, '.', 'Linewidth', 1, 'Color', [0.8 0.8 0.8], 'MarkerSize', 14);
xlabel({'Profondeur à partir'; 'de la surface de la peau  [mm]'}); ylabel('Surface activée  [mm^2]');
%legend('Fibres A\delta','Fibres C');
title('Surface de la peau activée');

subplot(2,2,4)
plot(x, dt*l, 'k', 'Linewidth', 1); hold on
plot(x, dt*m, 'Linewidth', 1, 'Color', [0.8 0.8 0.8]);
plot(x, dt*l, '.k', 'MarkerSize', 14); hold on
plot(x, dt*m, '.', 'Color', [0.8 0.8 0.8], 'MarkerSize', 14);
xlabel({'Profondeur à partir'; 'de la surface de la peau  [mm]'}); ylabel('Durée  [s]');
title('Durée d''activation');
legend('Fibres A\delta','Fibres C');