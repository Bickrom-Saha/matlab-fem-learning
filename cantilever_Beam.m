clear;
clc;
close all;

% parameters

L = 1; % Length in m
E = 200e9; %young's modulus in Pa
b = 0.02; % cross section width in m
h = 0.04; % cross section height in m
P = 100; % downward tip force in N

I = b*h^3/12; % second moment of area in m^4

% analytical check

w_tip_exact = -P*L^3/(3*I*E); % Euler Bernoulli beam
theta_tip_exact = -P*L^2/ (2*I*E);

fprintf('Exact tip displacement = %.6f mm\n',w_tip_exact*1000);
fprintf('Exact Tip rotation = %.6f radian\n',theta_tip_exact);

% Mesh: 4 2-node beam elements (Linear)
nElem = 4;
nNode = nElem +1; 
nodeCoords = linspace(0,L,nNode)'; % node coordinates
elemConn = [(1:nElem).' (2:nNode).']; % element connectivity
disp(elemConn);

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
    K(dofs,dofs) = K(dofs,dofs) + ke;
end

disp('Global stiffness matrix size: ');
disp(size(K));
disp(K);

% downward force at the final displacement DOF 
tipDOF = nDOF - 1; % and nDOF is the final rotational dof
F(tipDOF) = -P;

% Fixed left end (w1 =0, theta1 = 0)
fixedDofs = [1, 2];
alldofs = 1:nDOF;
freedofs = setdiff(alldofs,fixedDofs);
U = zeros(nDOF , 1);
U(freedofs) = K(freedofs,freedofs) \ F(freedofs);

%extract displacement and roation
w = U(1:2:end); % entries = 1, 3, 5...
theta = U(2:2:end);  % entries = 2, 4, 6...

% support reactions
R = K*U - F;
fprintf('FEM tip displacement = %.6f mm\n', w(end)*1000 );
fprintf('FEM tip rotation = %.6f rad\n', theta(end));
fprintf('FEM support reaction = %.6f N\n',R(1));
fprintf('FEM support moment = %.6f N*m\n',R(2));

fprintf('Exact support reaction = %.6f N\n',P);

% plot deflected shape
figure;
plot(nodeCoords,w*1000, LineStyle="-", Marker="o", LineWidth=1.5, MarkerSize=6);

% plot the interpolated FEM deflection

figure;
hold on;

for e =1:nElem
    nodes = elemConn(e,:);
    i = nodes(1);
    j = nodes(2);
    Le = nodeCoords(j) - nodeCoords(i);
    dofs = [2*i-1, 2*i, 2*j-1, 2*j];
    Ue = U(dofs);

    % sample positions within this element
    s = linspace(0,1,30).';

    % Cubic Hermite shape fuction
    N1 = 1 - 3*s.^2 + 2*s.^3;
    N2 = Le*(s - 2*s.^2 + s.^3);
    N3 = 3*s.^2 - 2*s.^3;
    N4 = Le*(-s.^2 + s.^3);

    N = [N1 N2 N3 N4];

    xplot = nodeCoords(i) +Le*s;  % s = (x - xi) / Le
    wplot = N*Ue;

    plot(xplot,wplot*1000, LineWidth=2, DisplayName= sprintf('Element %d', e));
end

% Analytical defelction at closely spaced positions
xExact = linspace(0,L,200);
wExact = -P*xExact.^2 .*(3*L - xExact)/(6*I*E);
plot(xExact,wExact*1000, 'r--', LineWidth=2, DisplayName="Analytical");
plot(nodeCoords, w*1000, 'bo', LineWidth= 2, DisplayName= 'FEM nodes', MarkerFaceColor='b');

xlabel('position along beam, m');
ylabel('Transverse Deflection, w(x) (mm)');
title('Cantilever Beam Deflection under a tip load');
legend(Location="best");
axis tight;
hold off;

% Recover Bending moment and bottom surface stress

figure;

for e =1:nElem
    nodes = elemConn(e,:);
    i = nodes(1);
    j = nodes(2);
    Le = nodeCoords(j) - nodeCoords(i);
    dofs = [2*i-1, 2*i, 2*j-1 , 2*j];
    Ue = U(dofs);

    s =linspace(0,1,30);
    xplot = nodeCoords(i) + Le*s;
   
    syms xi

    Nsym = [1 - 3*xi^2 + 2*xi^3, ...
            Le*(xi - 2*xi^2 + xi^3), ...
            3*xi^2 - 2*xi^3, ...
            Le*(-xi^2 + xi^3)];

    B_curvature = double(subs(diff(Nsym, xi, 2), xi, s.')) / Le^2;
   % disp(B_curvature);
    curvature = B_curvature*Ue;
    moment  = E*I*curvature;

    % sigma_xx = -M*y/I
    %Bottom surface : y = -h/2

    stress_bottom = moment*h/(2*I);
    subplot(2,1,1);
    hold on;
    plot(xplot, moment, '-', LineWidth=1.5, DisplayName= sprintf('FEM element %d',e));
    legend(Location='best');
    subplot(2,1,2);
    hold on;
    plot(xplot, stress_bottom/1e6, '-', LineWidth=1.5, DisplayName= sprintf('FEM element %d',e));
    legend(location='best');

end

% analytical compariseon
xExact = linspace(0,L,200);
moment_exact = -P*(L - xExact);
stress_bottom_exact = moment_exact*h/(2*I);

subplot(2,1,1);
plot(xExact, moment_exact, 'g--', LineWidth= 1.5, DisplayName='analytical');
xlabel('Position, m');
ylabel('Bending moment, N*m');
title('Bending moment along the beam');
legend(Location='best');
axis tight;
grid on;
hold off

subplot(2,1,2);
plot(xExact, stress_bottom_exact/1e6, 'g--', LineWidth= 1.5, DisplayName='analytical');
xlabel('Position, m');
ylabel('Bottom-surface stress, MPa');
title('Bottom surface axail stress');
grid on;
legend(Location='best');
axis tight;
hold off;