clear all;
close all;

nPOrder = 5;
nMOrder = 3;

p   = [1:nPOrder];
int = zeros(nPOrder, nMOrder);
sol = ones(nPOrder, 1) * -0.3926;

int(1, 1) = -0.2886;
int(2, 1) = -0.3071;
int(3, 1) = -0.3123;
int(4, 1) = -0.3130;
int(5, 1) = -0.3131;

int(1, 2) = -0.3760;
int(2, 2) = -0.3913;
int(3, 2) = -0.3921;
int(4, 2) = -0.3920;
int(5, 2) = -0.3920;

int(1, 3) = -0.3762;
int(2, 3) = -0.3922;
int(3, 3) = -0.3922;
int(4, 3) = -0.3928;
int(5, 3) = -0.3928;

plot(p, [int sol], '*-');