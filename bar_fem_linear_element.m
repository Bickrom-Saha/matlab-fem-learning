clear;
clc;
close all;
% problem inputs: m,N,Pa
L= 1; % bar length 1m
A= 100e-6; % cross sectional area in m^2
E= 200e9; %youngs modulus in Pa
P= 1000; % Axial force applied on the right end of the bar in N

% Analytical solution for checking out solution later
u_exact = P*L/(A*E);

fprintf('Exact end displacement = %.5f mm\n', u_exact*1000);

% Mesh: divide the bar into two equal elements.
nElem = 2;
nNode = nElem + 1;

% Lenght of each element in this uniform mesh: 
Le = L/nElem;

% Coordinates of the nodes:
x= linspace(0,L,nNode); % intial 0, final L , devides into n points. 

% Connectivity: the two node numbers belonging to each element
elements = [1 2;
            2 3];

disp('Node coordinates [m]: ');
disp(x);

disp('Element Connectivity: ');
disp(elements);

% Linear Element stiffness matrix: bar ( page 65 Dr. Clayton Pettit)
ke = (E*A/Le)*[1 -1;
               -1 1];
disp('Element stiffness matrix: ');
disp(ke);

% Initialize the global stiffness matrix:
K = zeros(nNode, nNode);

% Assemble each element into the global matrix:
for e = 1:nElem
    nodes = elements(e, :);

    K(nodes, nodes) = K(nodes, nodes) + ke;
end
disp('Global Stiffness Matrix: ');
disp(K);

% Applied nodal forces: 
F = zeros(nNode, 1);
F(3) = P;

% Prescribed displacements and unknown displacements: 
u = zeros(nNode, 1);
fixed = 1;
free = setdiff(1:nNode, fixed); % setdiff([1 2 3] , 1) gives free = [2 3]

%solve for the unknown displacements:
u(free) = K(free, free) \ F(free); % x = A\B means find x such that A*x = B
                                   % u([2 3]) = K([2 3],[2 3]) \ F([2 3])
                                   % u(rows) = K(rows, colums) \ F(rows)

disp('Nodal displacemets : ');
disp(u*1000);

% Support reaction and residual forces: 

R = (K * u) - F;    % page 67 dr. clayton
disp('Reaction and residual forces: ');
disp(R);

% Calculate element stress and strain:
strain = zeros(nElem, 1);
stress = zeros(nElem, 1);

for e = 1:nElem
    nodes = elements(e, :);
    extension = u(nodes(2)) - u(nodes(1));
    strain(e) = extension/ Le;
    stress(e) = E*strain(e);
end
disp('Element strains: ');
disp(strain);
disp('element stresses [MPa]: ');
disp(stress/ 1e6);

figure;
plot(x, u*1000, '-o', 'LineWidth',1.5);

xlabel('Position along bar, m');
ylabel('Axial displacement, mm');
title('Axial displacement along the length of the bar');
grid on; 

% Scale from zero load to the current solved load
loadfactor = linspace(0,1,21);

% Linear elastic loading response
straincurve = loadfactor * strain(1);
stresscurve = loadfactor * stress(1);
tipdisplacement = loadfactor * u(end);
appliedforce = loadfactor * P;

% stress vs strain
figure;
plot(straincurve, stresscurve);
xlabel('strain');
ylabel('stress');
grid on;

% force vs displacement
figure;
plot(tipdisplacement, appliedforce);
xlabel('tip displacement');
ylabel('applied force');
grid on;