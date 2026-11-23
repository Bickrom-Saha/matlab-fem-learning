% Simply supported beam under uniform load

clc;
clear;
close all;

% parameters

L = 1; % Length in m
E = 200e9; %young's modulus in Pa
b = 0.02; % cross section width in m
h = 0.04; % cross section height in m
q0 = 100; % downward load magnitude, N/m
qy = -q0; % signed transverse load, N/m
I = b*h^3/12; % second moment of area in m^4

% Mesh: 4 2-node beam elements (Linear)
nElem = 40;
nNode = nElem +1; 
nodeCoords = linspace(0,L,nNode)'; % node coordinates
elemConn = [(1:nElem).' (2:nNode).']; % element connectivity
%disp(elemConn);

% 2 DOF per node
% Global ordering [ w1, theta1, w2, theta2....]
nDOF = 2*nNode;

K = zeros(nDOF);
F = zeros(nDOF,1);

% element stiffness and assembly
for e = 1:nElem
    nodes = elemConn(e,:);

    i = nodes(1);
    j = nodes(2);

    Le = nodeCoords(j) - nodeCoords(i);
    ke = (E*I/Le^3) *[12,   6*Le,   -12,    6*Le;
        6*Le, 4*Le^2, -6*Le,  2*Le^2;
        -12,  -6*Le,  12,     -6*Le;
        6*Le, 2*Le^2, -6*Le,  4*Le^2];

    % Global DOFs ; [w_i, theta_i, w_j, theta_j]
    dofs = [2*i-1, 2*i, 2*j-1, 2*j];
    fe = qy * [Le/2;
        Le^2/12;
        Le/2;
        -Le^2/12];
    K(dofs,dofs) = K(dofs,dofs) + ke;
    F(dofs) = F(dofs) + fe;
end

disp('Global stiffness matrix size: ');
disp(size(K));
%disp(K);

% Fixed Transverse displacement at both end
% Leave both end rotation free

fixed = [1, nDOF-1];
free = setdiff(1:nDOF, fixed);

% Solve
U = zeros(nDOF,1);
U(free) = K(free,free)\F(free);

w = U(1:2:end); % return 1, 3, 5...
theta = U(2:2:end); % returns 2, 4, 6...

R = K*U -F;

fprintf('Left support reaction = %.6f N\n',R(1));
fprintf('Right support reaction = %.6f N\n', R(end-1));

xExact = linspace(0,L,300);
wExact= -q0*xExact.*(L^3 - 2*L*xExact.^2 + xExact.^3)/(24*E*I);

MExact = q0*xExact.*(L-xExact)/2;

figure;
subplot(2,1,1);
plot(xExact, wExact*1000, 'r-', LineWidth=1.5, DisplayName='Analytical');
hold on;
subplot(2,1,2);
plot(xExact, MExact, 'g--', LineWidth=1.5, DisplayName='analytical');
hold on;

%Evaluate FEM results inside each elements

for e = 1:nElem
    nodes = elemConn(e, :);
    i = nodes(1);
    j = nodes(2);
    Le = nodeCoords(j) - nodeCoords(i);

    dofs = [2*i-1, 2*i, 2*j-1, 2*j];
    Ue = U(dofs);
    s =linspace(0,1,50).';
    xplot = nodeCoords(i) +Le*s;
    % Displacement interpolation 
    N1 = 1 - 3*s.^2 + 2*s.^3;
    N2 = Le*(s - 2*s.^2 + s.^3);
    N3 = 3*s.^2 - 2*s.^3;
    N4 = Le*(-s.^2 + s.^3);

    N = [N1 N2 N3 N4];

    syms xi

    Nsym = [1 - 3*xi^2 + 2*xi^3, ...
        Le*(xi - 2*xi^2 + xi^3), ...
        3*xi^2 - 2*xi^3, ...
        Le*(-xi^2 + xi^3)];

    B_curvature = double(subs(diff(Nsym, xi, 2), xi, s)) / Le^2;
    wFEM = N*Ue;
    MFEM = E*I*(B_curvature*Ue);

    if e == nElem/2
        Bmidspan = double(subs(diff(Nsym, xi, 2), xi, 1)) / Le^2;
        Mmidspan = E*I*(Bmidspan*Ue);
    end

    % Label only the first FEM segment in the legend

    visibility ='off';
    if e == 1
        visibility ='on';
    end

    subplot(2,1,1);
    plot(xplot, wFEM*1000, 'b:', LineWidth=1.5, DisplayName='FEM', HandleVisibility=visibility);

    subplot(2,1,2);
    plot(xplot, MFEM, 'b:', LineWidth=1.5, DisplayName='FEM', HandleVisibility=visibility);

end

subplot(2,1,1);
legend(Location='best');
xlabel('Distance along beam, x (m)');
ylabel('Deflection, w(x) (mm)');
title('Deflection under distributed loading');
grid on;
hold off;

subplot(2,1,2);
legend(Location='best');
xlabel('Distance along beam, x (m)');
ylabel('Bending moment, M(x) (N.m)');
title('Bending moment under distributed loading');
grid on; 
hold off;

midNode = nElem/2 +1; % number of Element is even
fprintf('Midspan deflection (FEM) = %.6f mm\n',w(midNode)*1000);
fprintf('Midspan moment (FEM) = %.6f Nm\n',Mmidspan);