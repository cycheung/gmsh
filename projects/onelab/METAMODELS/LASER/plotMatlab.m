clf
clf('reset')

%temperatures at x=0 as a function of (time,depth)     
A = importdata('temp.txt');

Tad = 320; %temp activation fivres Ad
nbz = 5;
z = [0 0.05 0.10 0.15 0.20];

z0=2;
z1=8;
z2=14;
z3=20;
z4=26;

%maximum temp and index I of when it is reached I = 1 .. nbTimeStep
I = ones(5,1);
[m0,I(1)]= max(A(:,z0));
[m1,I(2)]= max(A(:,z1));
[m2,I(3)]= max(A(:,z2));
[m3,I(4)]= max(A(:,z3));
[m4,I(5)]= max(A(:,z4));

I 

%time step
disp('size tempx=0, dt, nbTime')
size(A)
dt = A(2,1)-A(1,1)
nbTimeStep = size(A,1)

%temperature at z= 0 as a function of x
B = importdata('tempsurf.txt');

%temperatures at differents depths as a function of time
B0 = importdata('temp0.txt');
B1 = importdata('temp1.txt');
B2 = importdata('temp2.txt');
B3 = importdata('temp3.txt');
B4 = importdata('temp4.txt');

BBp = [B0', B1', B2', B3', B4'];
BB = BBp';

  disp('sizes BBBB')
size(BB)
size(B0)

nbRow = size(B0,1)
nbCol = size(B0,2)
nbDT = nbRow/nbTimeStep
for i=1:nbz
  ind = I(i);
  disp('sizes')
  1+nbDT*(ind-1)
  1+nbDT*(ind)
  (i-1)*nbCol+5
  size(BB)
  BB_x = BB(1+nbDT*(ind-1):nbDT*ind,(i-1)*nbCol+5);
  BB_T = BB(1+nbDT*(ind-1):nbDT*ind,(i-1)*nbCol+8);
  size(BB_T)
  Xhot = 0
  for j =1:nbDT
    if ( BB_T(j) > Tad)
       Xhot = BB_x(j);
    end
  end
  Aactive(i) = (pi*Xhot^2)*10^6
end

%plot
subplot(2,2,[1 2])
plot( A(:,1), A(:,z0), 'k', 'Linewidth', 1.5); hold on
plot( A(:,1), A(:,z1), 'k', 'Linewidth', 1.5); hold on
plot( A(:,1), A(:,z2), 'k', 'Linewidth', 1.5); hold on
plot( A(:,1), A(:,z3), 'k', 'Linewidth', 1.5); hold on
plot( A(:,1), A(:,z4), 'k', 'Linewidth', 1.5); hold on
plot( A(:,1), 320*ones(size(A(:,1))), 'LineStyle', '-', 'Color', [0.8 0.8 0.8], 'Linewidth', 1); hold on

legend('z = 0 mm', 'z = 0.05 mm', 'z = 0.10 mm', 'z = 0.15 mm', 'z = 0.20 mm', 'Seuil fibres A\delta', 'Seuil fibres C');
xlabel('Temps  [s]'); ylabel('TempÃ©rature  [K]');
title('Température de la peau exposée à  un laser CO2 pour différentes profondeurs');

subplot(2,2,3)
plot( B(:,5), B(:,8), 'k', 'Linewidth', 1.5); hold on
xlabel('x  [s]'); ylabel('Température  [K]');
title('Température à la surface de la peau');


subplot(2,2,4)
plot( z, Aactive, '.k', 'Linewidth', 1, 'MarkerSize', 14); hold on
xlabel({'Profondeur ‡ partir'; 'de la surface de la peau  [mm]'}); ylabel('Surface activÈe  [mm^2]');
title('Surface de la peau activÈe');

