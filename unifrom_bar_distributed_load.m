clear;
clc;
close all;

% Problem statement: solve the unifrom bar with a distributed axial load,
% using quadratic elements and gaussian quadrature. left end is fixed

%parameters
L = 1;      % bar length m
E = 200e9;  %youngs modulus Pa
A = 100e-6; % cross sectional area m^2
q = 1000;   % Axial distributed load N/m

% Mesh: 2 elements
nElem = 2;
nNode = 2*nElem +1;
x = linspace(0,L,nNode);

% local node order: 
elements = [1 2 3;
            3 4 5];
% Two-point Gaussian Quadrature: 2 degree polynomial
gaussWeights = [1 1];
gaussPoints = [-1/sqrt(3) 1/sqrt(3)];

% Initialize global matrix
K = zeros(nNode,nNode);
F = zeros(nNode,1);

% Element integration and assembly
for e = 1:nElem
    nodes = elements(e, :);

    % Element coordinates as a column vector
    xe = x(nodes).';
    ke = zeros(3, 3);
    fe = zeros(3,1);

    for g = 1: length(gaussPoints)
        xi = gaussPoints(g);
        weight = gaussWeights(g);
        % quadratic shape function
        N = [0.5*xi*(xi-1), 1-xi^2, 0.5*xi*(xi+1)];

        % derivatives with respect to natural coordinate xi
        dN_dxi = [(xi-0.5), -2*xi, (xi+0.5)];

        % Jacobian dx/dxi
        J = dN_dxi * xe;

        % Derivatives with respect to physical coordinate x
        B = dN_dxi / J;
        
        % Add this gauss point's contribution
        ke = ke + (B.' * E * A * B) * J* weight;
        fe = fe + (N.' * q) * J * weight;
    end

    % Assemble into the global system
    K(nodes, nodes) = K(nodes, nodes) + ke;
    F(nodes) = F(nodes) + fe;
end

% Boundary conditions and sulution
u = zeros(nNode, 1);
fixed = 1;
free = setdiff(1:nNode, fixed);
u(free) = K(free, free) \ F(free);
R = K*u - F

% Analytical Solution at the nodes
u_exact = (q/(E*A))* (L*x - x.^2/2);

disp('Applied equivalent Nodal loads,(N) : ');
disp(F.');

disp('FEM Computed displacements, u (mm): ');
disp(u.');

disp('Analytical displacements, u_exact: ');
disp(u_exact.'*1000);

fprintf('Total applied load = %0.3f N\n', sum(F));
fprintf('Support reaction = %0.3f N\n', R(1));  

figure;
plot(x, F, '-o', 'LineWidth', 1.5);
xlabel('Position (m)');
ylabel('Nodal load (N)');
title('Position vs. Nodal Load');
grid on;
figure;
plot(x, u * 1000, 'bo-', LineWidth= 1.5, DisplayName= sprintf('FEM'));
hold on;
plot(x, u_exact * 1000, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Analytical');
xlabel('Position (m)');
ylabel('Displacement (mm)');
title('FEM vs. Analytical Displacement');
legend('Location', 'best');
grid on;
hold off;
% Exact solution using dsolve
syms  u(x); % x no became symbol not variable.
syms E A L q positive; %same here, symbols
% governing equation
ODE = E*A*diff(u,x,2)+q == 0;
du = diff(u,x);
BC = [u(0) == 0, du(L) == 0];
uSol = simplify(dsolve(ODE, BC));
disp(uSol);
x_nodes = linspace(0, 1, nNode);
uSol_expr = formula(uSol);
uSol_nodes = subs(uSol_expr, x, x_nodes);
uSol_nodes = subs(uSol_nodes, {L, A, E, q}, {1, 100e-6, 200e9, 1000});
uSol_nodes = double(uSol_nodes);

disp('Analytical solution at each node (mm):');
disp(uSol_nodes.'*1000);
figure;
plot(x_nodes, uSol_nodes * 1000, 'ro', 'LineWidth', 1.5, 'DisplayName', 'Analytical');
hold on;
