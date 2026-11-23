% Quadratic  elements have 3 nodes for each element
% different cross sectional area
clear;
clc;
close all;
% problem inputs: m,N,Pa
L= 1; % bar length 1m
A= [100e-6, 50e-6]; % cross sectional area in m^2
E= 200e9; %youngs modulus in Pa
P= 2000; % Axial force applied on the right end of the bar in N

% Analytical solution for checking out solution later
u_exact = P*(L/2)/(A(1)*E) + P*(L/2)/(A(2)*E);

fprintf('Exact end displacement = %.5f mm\n', u_exact*1000);

% Mesh: divide the bar into two equal elements.
nElem = 2;
nNode = 2*nElem + 1;

% Coordinates of the nodes:
x= linspace(0,L,nNode); % intial 0, final L , devides into n points. 

% Connectivity: the two node numbers belonging to each element
elements = [1 2 3;
            3 4 5];

disp('Node coordinates [m]: ');
disp(x);

disp('Element Connectivity: ');
disp(elements);

% Initialize the global stiffness matrix:
K = zeros(nNode, nNode);

% Assemble each element into the global matrix:
for e = 1:nElem
    nodes = elements(e, :);
    % element length and stiffness
    Le = x(nodes(3))- x(nodes(1));
    ke = (E*A(e)/(3*Le))*[7 -8 1;
                       -8 16 -8;
                       1  -8 7];

    K(nodes, nodes) = K(nodes, nodes) + ke;
end

disp('Element stiffness matrix: ');
disp(ke);
disp('Global Stiffness Matrix: ');
disp(K);

% Applied nodal forces: 
F = zeros(nNode, 1);
F(end) = P;

% Prescribed displacements and unknown displacements: 
u = zeros(nNode, 1);
fixed = 1;
free = setdiff(1:nNode, fixed); % setdiff([1 2 3] , 1) gives free = [2 3]

%solve for the unknown displacements:
u(free) = K(free, free) \ F(free); % x = A\B means find x such that A*x = B
                                   % u([2 3]) = K([2 3],[2 3]) \ F([2 3])
                                   % u(rows) = K(rows, colums) \ F(rows)

disp('Nodal displacemets [mm] : ');
disp(u*1000);

% Support reaction and residual forces: 

R = (K * u) - F;    % page 67 dr. clayton
disp('Reaction and residual forces [N]: ');
disp(R);

% Calculate element stress and strain:
strain = zeros(nElem, 1);
stress = zeros(nElem, 1);

for e = 1:nElem
    nodes = elements(e, :);
    Le = x(nodes(3))-x(nodes(1));
    Bmid = [-1 0 1]/Le;
    strain(e) = Bmid * u(nodes);
    stress(e) = E*strain(e);
end
disp('Element midpoint strains: ');
disp(strain);
disp('element midpoint stresses [MPa]: ');
disp(stress/ 1e6);

figure;
plot(x, u*1000, '-o', 'LineWidth',1.5);

xlabel('Position along bar, m');
ylabel('Axial displacement, mm');
title('Axial displacement along the length of the bar');
grid on; 

% Natural coordinate of the left , middle and right nodes
xiNodes = [-1,0,1];

% one row per element, one column per local node
nodalStrain = zeros(nElem, 3);
nodalStress = zeros(nElem, 3);

for e = 1:nElem
    nodes = elements(e, :);
    Le = x(nodes(3)) - x(nodes(1));
    ue = u(nodes);
    for a =1:3
        xi = xiNodes(a);
        % Shape- Function derivatives with respect to x
        B = (2 / Le) * [xi - 0.5, -2*xi, xi+0.5];

        nodalStrain(e, a) = B* ue;
        nodalStress(e, a) = E * nodalStrain(e, a);
    end
end

disp('strains at each element''s thress Nodes: ');
disp(nodalStrain);
disp(' stresses at each elemet''s three nodes [MPa]: ');
disp(nodalStress/ 1e6);

% plot each element separately to preserve the interface jump 
figure;
subplot(2,1,1);
hold on;
for e = 1:nElem
    nodes = elements(e, :);
    plot(x(nodes), nodalStrain(e, :), '-o', LineWidth= 1.5, DisplayName= sprintf('Element %d',e));
end
xlabel('Position [m]');
ylabel('Axial Strain [-]');
title("strain at element nodes");
legend('show');
grid on;
hold off;

subplot(2,1,2);
hold on;
for e=1:nElem
    nodes = elements(e, :);
    plot(x(nodes), nodalStress(e, :)/1e6, '-o', LineWidth= 1.5, DisplayName= sprintf('Element %d',e));
end
xlabel('Position [m]');
ylabel('Axial Stress [MPa]');
title("stress at element nodes");
legend('show');
grid on;
hold off;
