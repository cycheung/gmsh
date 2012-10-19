close all;
clear all;

%% l2 [Order][Mesh]

%% Sin(10x) + Sin(10y) + Sin(10z)

l2 = ...
    [...
        +1.229975e+00 , +7.127187e-01 , +2.574109e-01 , +6.358030e-02 , +1.398086e-02 ; ...
        +1.121877e+00 , +2.585136e-01 , +5.192648e-02 , +8.839964e-03 , +1.373343e-03 ; ...
    ];


h = [1, 1/2, 1/4, 1/8, 1/16];
p = [1:2];

P = size(p, 2);
H = size(h, 2);

delta = zeros(P, H - 1);

for i = 1:H-1
    delta(:, i) = ...
        (log10(l2(:, i + 1)) - log10(l2(:, i))) / ...
        (log10(1/h(i + 1))   - log10(1/h(i)));
end

delta

figure;
loglog(1./h, l2, '-*');
grid;
