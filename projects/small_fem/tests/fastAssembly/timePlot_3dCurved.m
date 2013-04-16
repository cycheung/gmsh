%% Clear
close all;
clear all;

%% Get Data
data = dlmread("timing_3dCurved.csv", ";");

order    = data(1, 2:end);
unknowns = data(2, 2:end);
slow     = data(3, 2:end);
fast     = data(4, 2:end);

%% Time
figure;
hold on;
semilogy(order, slow, '-or');%, 0.125);
semilogy(order, fast, '-ob');%, 0.125);
hold off;

title("Fast and Slow Assembly Algorithm");
xlabel("Interplation Order [-]");
ylabel("Time [s]");

legend([{"Slow"},
        {"Fast"}], 'location', 'northwest');
grid;

%print 'time.svg'

%% Unknowns
figure;
plot(order, unknowns);
title("Number of Unknowns");
xlabel("Interplation Order [-]");
ylabel("Number of Unknowns [-]");

legend([{"Number of Unknowns"}]);
grid;
