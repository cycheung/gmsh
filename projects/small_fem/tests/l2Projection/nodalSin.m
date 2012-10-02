close all;
clear all;

%% l2 [Order][Power]

%% Sin(10x) + Sin(10y)
l2 = ...
    [...
        +9.824759e-01 , +1.228302e+00 , +9.391098e-01 , +3.471288e-01 , +6.947808e-02 ; ...
        +1.049055e+00 , +1.103665e+00 , +3.087024e-01 , +5.534351e-02 , +9.315162e-03 ; ...
        +1.162833e+00 , +9.113153e-01 , +1.194850e-01 , +8.849612e-03 , +5.227025e-04 ; ...
        +1.137504e+00 , +3.108951e-01 , +2.401343e-02 , +1.037393e-03 , +3.936262e-05 ; ...
        +1.269043e+00 , +2.610065e-01 , +6.347677e-03 , +1.054095e-04 , +1.652859e-06 ; ...
        +9.831276e-01 , +5.001032e-02 , +8.915199e-04 , +9.697960e-06 , +8.755089e-08 ; ...
        +9.255320e-01 , +3.346505e-02 , +1.795654e-04 , +7.042222e-07 , +2.808220e-09 ; ...
        +5.770821e-01 , +3.946435e-03 , +1.897450e-05 , +5.129431e-08 , +1.145673e-10 ; ...
    ];

%% Sin(2x) + Sin(2y)
%l2 = ...
%    [...
%        +3.023927e-01 , +2.019098e-01 , +4.088128e-02 , +9.640896e-03 , +2.368553e-03 ; ...
%        +1.848512e-01 , +2.225958e-02 , +3.927832e-03 , +5.514234e-04 , +7.203640e-05 ; ...
%        +1.488307e-02 , +3.737334e-03 , +2.051181e-04 , +1.245317e-05 , +7.656078e-07 ; ...
%        +9.619936e-03 , +2.341946e-04 , +9.749372e-06 , +3.320698e-07 , +1.073575e-08 ; ...
%        +3.209790e-04 , +2.975511e-05 , +4.232692e-07 , +6.529912e-09 , +9.973727e-11 ; ...
%        +2.386183e-04 , +1.311533e-06 , +1.308600e-08 , +1.102963e-10 , +1.866212e-12 ; ...
%        +5.381028e-06 , +1.295047e-07 , +4.691249e-10 , +2.854052e-12 , +3.382652e-12 ; ...
%        +3.224924e-06 , +4.302990e-09 , +1.108367e-11 , +6.350352e-12 , +1.250808e-11 ; ...
%    ];

m = [2, 8, 32, 128, 512];
o = [1:8];
O = size(o, 2);
M = size(m, 2);

delta = zeros(O, M - 1);

for i = 1:M-1
    delta(:, i) = ...
        (log10(l2(:, i + 1)) - log10(l2(:, i))) / ...
        (log10(m(i + 1))     - log10(m(i)));
end

delta

figure;
loglog(m, l2, '-*');
grid;