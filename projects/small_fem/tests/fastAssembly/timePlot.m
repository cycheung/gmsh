%% Clear
close all;
clear all;

%% Get Data
data = dlmread("timing.csv", ";");

% Order
order    = data(1, 3:end);
unknowns = data(2, 3:end);

% Slow
sPrecmpt = data(3, 3:end);
sAlloc   = data(4, 3:end);
sDof     = data(5, 3:end);
sTerm    = data(6, 3:end);
sAdd     = data(7, 3:end);
slow     = [sPrecmpt' sAlloc' sDof' sTerm' sAdd'];

% Fast
fPrecmpt = data(8,  3:end);
fAlloc   = data(9,  3:end);
fDof     = data(10, 3:end);
fTerm    = data(11, 3:end);
fAdd     = data(12, 3:end);
fast     = [fPrecmpt' fAlloc' fDof' fTerm' fAdd'];

%% Time
figure;
hold on;
bar(order - 0.125, slow, 0.125, 'stacked');
bar(order + 0.125, fast, 0.125, 'stacked');
hold off;

title("Fast and Slow Assembly Algorithm");
xlabel("Interplation Order [-]");
ylabel("Time [s]");

legend([{"Precomputing"},
        {"Preallocation"},
        {"Dof Lookup"},
        {"Term Computation"},
        {"Term Addition"}], 'location', 'northwest');
grid;

%% Unknowns
figure;
plot(order, unknowns);
title("Number of Unknowns");
xlabel("Interplation Order [-]");
ylabel("Number of Unknowns [-]");

legend([{"Number of Unknowns"}]);
grid;
