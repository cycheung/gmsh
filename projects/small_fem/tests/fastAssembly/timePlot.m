%% Clear
close all;
clear all;

%% Get Data
data = dlmread("logTimeBlas.csv", ";");

off = 0; % Stops before end of data

order    = data(2:end - off, 1);
unknowns = data(2:end - off, 2);
fast     = data(2:end - off, 3);
slow     = data(2:end - off, 4);

%% Time
figure;
bar(order, [fast slow]);
title("Fast and Slow Assembly Algorithm");
xlabel("Interplation Order [-]");
ylabel("Time [s]");

legend([{"With BLAS 3 Operations"},
        {"Without BLAS 3 Operations"}]);
grid;

%% Unknowns
figure;
plot(order, unknowns);
title("Number of Unknowns");
xlabel("Interplation Order [-]");
ylabel("Number of Unknowns [-]");

legend([{"Number of Unknowns"}]);
grid;
