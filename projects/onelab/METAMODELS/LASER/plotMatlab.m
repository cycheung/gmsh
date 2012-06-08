clf
clf('reset')

%temperatures at x=0 as a function of (time,depth)     
S = importdata('temp.txt');

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
[m0,I(1)]= max(S(:,z0));
[m1,I(2)]= max(S(:,z1));
[m2,I(3)]= max(S(:,z2));
[m3,I(4)]= max(S(:,z3));
[m4,I(5)]= max(S(:,z4));

%time step
disp('dt, nbTime')
dt = S(2,1)-S(1,1)
nbTimeStep = size(S,1)

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

size(BB)

nbRow = size(B0,1);
nbCol = size(B0,2);
nbDT = nbRow/nbTimeStep;
for i=1:nbz
  ind = I(i);
  ind_bb = nbRow*(i-1);
  BB_x = BB(ind_bb+nbDT*(ind-1)+1:ind_bb+nbDT*ind,5);
  BB_T = BB(ind_bb+nbDT*(ind-1)+1:ind_bb+nbDT*ind,8);
  Xhot = 0;
  for j =1:nbDT
    if ( BB_T(j) > Tad)
       Xhot = BB_x(j);
    end
  end
  Aactive(i) = (pi*Xhot^2)*10^6; %active surface in mm
end

Aactive

%active surface for different time
A0 = importdata('active0.txt');
A1 = importdata('active1.txt');
A2 = importdata('active2.txt');
A3 = importdata('active3.txt');
A4 = importdata('active4.txt');
ActiveSurf = [max(A0(:,8)), max(A1(:,8)),max(A2(:,8)),max(A3(:,8)),max(A4(:,8))]*10^6;

ActiveSurf 

%plot
subplot(3,2,[1 2])
plot( S(:,1), S(:,z0), 'k', 'Linewidth', 1.5); hold on
plot( S(:,1), S(:,z1), 'b', 'Linewidth', 1.5); hold on
plot( S(:,1), S(:,z2), 'r', 'Linewidth', 1.5); hold on
plot( S(:,1), S(:,z3), 'g', 'Linewidth', 1.5); hold on
plot( S(:,1), S(:,z4), 'y', 'Linewidth', 1.5); hold on
plot( S(:,1), 320*ones(size(S(:,1))), 'LineStyle', '-', 'Color', ...
      [0.8 0.8 0.8], 'Linewidth', 1); hold on
plot( S(:,1), 320*ones(size(S(:,1))), 'LineStyle', '-', 'Color', ...
      [0.8 0.8 0.8], 'Linewidth', 1); hold on
plot( S(:,1), 314*ones(size(S(:,1))), 'LineStyle', '-', 'Color', [0.8 0.8 0.8], 'Linewidth', 1); hold on
hLegend = legend('z0', 'z1', 'z2', 'z3', 'z4', 'Threshold fibre A\delta', 'Threshold fibre C');
hXLabel1 = xlabel('Time  [s]'); hYLabel1 = ylabel('Temperature  [K]');
hTitle1 = title ('Skin Temperature of a laser CO2 for different depths');

subplot(3,2,[3 4])
B0l = B0(nbDT*9+1:nbDT*10,:);
B1l = B1(nbDT*9+1:nbDT*10,:);
B2l = B2(nbDT*9+1:nbDT*10,:);
B3l = B3(nbDT*9+1:nbDT*10,:);
B4l = B4(nbDT*9+1:nbDT*10,:);
plot( B0l(:,5)*1000, B0l(:,8), 'k', 'Linewidth', 1.5); hold on
plot( B1l(:,5)*1000, B1l(:,8), 'b', 'Linewidth', 1.5); hold on
plot( B2l(:,5)*1000, B2l(:,8), 'r', 'Linewidth', 1.5); hold on
plot( B3l(:,5)*1000, B3l(:,8), 'g', 'Linewidth', 1.5); hold on
plot( B4l(:,5)*1000, B4l(:,8), 'y', 'Linewidth', 1.5); hold on
hXLabel2 = xlabel('x [mm]'); hYLabel2 = ylabel('Temperature  [K]');
hTitle2 = title('Temperature at different depths at Time=Tlaser');

subplot(3,2,5)
plot( z, ActiveSurf, '--rs', 'Linewidth', 1.5); hold on
hTitle3 = title('Active surface (integrate)');
hXLabel3 =xlabel({'Depth  [mm]'}); hYLabel3 =ylabel('Surface [mm^2]');

subplot(3,2,6)
plot( z, Aactive,  '--rs', 'LineWidth', 1.5);  hold on;
hTitle4 = title('Active surface (pi*r2)');
hXLabel4 =xlabel({'Depth  [mm]'}); hYLabel4 =ylabel('Surface [mm^2]');

set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([hTitle1, hTitle2, hTitle3,hTitle4, hXLabel1, hYLabel1, hXLabel2, hYLabel2,hXLabel3, hYLabel3,hXLabel4, hYLabel4], ...
    'FontName'   , 'AvantGarde');
set([hLegend, gca]             , ...
    'FontSize'   , 11           );
set([hXLabel1, hYLabel1, hXLabel2, hYLabel2,hXLabel3, hYLabel3,hXLabel4, hYLabel4]  , ...
    'FontSize'   , 12          );
set( [hTitle1, hTitle2, hTitle3,hTitle4] , ...
    'FontSize'   , 14          , ...
    'FontWeight' , 'bold'      );

set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'on'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'YTick'       , 0:500:2500, ...
  'LineWidth'   , 1         );
