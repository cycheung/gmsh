close all;
clear all;

%% Vertex
vN = 9; 
V  = 1:vN;

%% Vertex ID
vId = V;

%% Edge
E  = [[1 5]; [5 9]; [9 1]; 
      [1 8]; [8 9];
      [9 7]; [7 8];
      [7 4]; [4 8];
      [5 2]; [2 6]; [6 5];
      [9 6];
      [6 3]; [3 9];
      [7 3]];

eN = size(E, 1);

%% Edge ID
eId = zeros(1, eN);

for i = 1:eN
    eId(i) = E(i, 1) + E(i, 2) * vN;
end

%% Cell
cN = 8;
C  = 1:cN;

%% Cell ID
cId = zeros(1, cN);

for i = 1:cN
    cId(i) = C(i) * vN * vN;
end

%% Plot Numbering
figure;
plot(vId, ones(vN) * 1.5, 'xr', ...
     eId, ones(eN) * 2.0, 'xb', ...
     cId, ones(cN) * 2.5, 'xk');

id = [vId eId cId];

axis([min(id) - 1, max(id) + 1, ...
      1          , 3]);

legend('Vertex', 'Edge', 'Cell', ...
       'Location', 'NorthWest');
