% 2D Quadrilateral 4 node solid element. thin rectangular plate
% Horizontal displacement on the right edge and vertical displacement at the bottom left corner
% The rest of the left edge can move vertically, allowing poisson contraction (Plane stress)

close all;
clear;
clc;

% Geometry, Material, and Traction

L = 1; % Length, m
H = 0.2; % Height, m
t = 0.01; % Thickness, m

E = 200e9; % Young's modulus, Pa
nu = 0.3; % Poisson's ratio
T = 1e6; % Right edge Traction, Pa

% Nodes of the plate
coords = [0 0;  % node 1 : bottom left
          L 0;  % node 2: bottom right
          L H;  % node 3: top right
          0 H]; % node 4: top left

% Counter-clockwise element connectivity

elements = [1 2 3 4];

% Degree of freedom
% ordering: [u1, v1, u2, v2, u3, v3, u4, v4]

nNode = size(coords, 1);
nDof = 2 * nNode;

K = zeros(nDof, nDof);
F = zeros(nDof, 1);

% Plane stress constitutive matrix

D = E/(1-nu^2)* [1 nu 0;
                 nu 1 0;
                 0 0 (1-nu)/2];

% Analytical Checks

uRe = T*L/E; % u_Right_Exact
vTe = -nu*T*H/E; % v_Top_Exact
%uLe = -T*H/E; % u_Left_Exact
% vBe = 0; %v_Bottom_Exact

tF = T*H*t; % total force

fprintf('Right edge horicontal displacement = %0.6f mm\n', uRe*1000);
fprintf('Top edge vertical displacement = %0.6f mm\n', vTe*1000);
fprintf('Total applied horizontal force = %0.6f N\n', tF);

% Four Gaussian points: Two in each direction
gaussian_points = [-1 1]/sqrt(3);
gaussian_weights = [1 1];

ke = zeros(8 , 8);

% Coordinates of this element's nodes: 4X2 matrix

xe = coords(elements(1,:), :);
%disp(xe);
%disp(coords);

for g1 = 1:2
    for g2 = 1:2
        xi = gaussian_points(g1);
        eta = gaussian_points(g2);

        weight = gaussian_weights(g1)*gaussian_weights(g2);
% N = 0.25[(1-xi)(1-eta), (1+xi)(1-eta), (1+xi)(1+eta), (1-xi)(1+eta)
        % Row 1: Derivatives w.r.t xi
        % Row 2: Derivatives w.r.t eta
        dN_natural = 0.25* [-(1-eta), (1-eta), (1+eta), -(1+eta);
                          -(1-xi),  -(1+xi),  (1+xi), (1-xi)];
        
        % Jacobian Convention:
        % J = [dx/dxi   dy/dxi;
        %      dx/deta  dy/deta];

        J = dN_natural * xe;
        detJ = det(J);
        
        if detJ <= 0
            error("Jacobian determinant is non-positive");
        end

        % Row1: dN/dx
        % Row2: dN/dy

        dN_physical = J \ dN_natural;

        % Strain displacement Matrix

        B = zeros(3, 8);

        for a = 1:4
            dNdx = dN_physical(1,a);
            dNdy = dN_physical(2,a);
            B(:, 2*a-1:2*a) = [dNdx 0;
                               0 dNdy;
                               dNdy dNdx];
        end
        ke = ke + B.' * D * B * t * detJ * weight;
    end
end
% this model contains only one element, with DOFs 1 through 8

K = K + ke;

disp("Element stiffness matrix size: ");
disp(size(ke));

% right edge traction: node 2 and 3
% DOFs: [u1 v1 u2 v2 u3 v3 u4 v4]

F(3) = T*H*t/2;
F(5) = T*H*t/2;

% Suppots
% u1 = 0, v1 = 0, u4 = 0

fixedDofs = [1 2 7];
alldofs = 1:nDof;
freedofs = setdiff(alldofs, fixedDofs);

% Solve

U = zeros(nDof,1);
U(freedofs) = K(freedofs, freedofs) \ F(freedofs);

% Extract displacement components:

u = U(1:2:end);
v = U(2:2:end);

% Reactions at Constrained DOFs and Residuals elsewhere:

R = K * U - F;

% Display results:

disp('  Node       u(mm)       v(mm)');
disp([(1:nNode).', u*1000, v*1000]);

fprintf('Total Horizontal Reaction = %.6f N\n', sum(R(1:2:end)));
fprintf('Total Vertical Reaction = %.6f N\n', sum(R(2:2:end)));

% Now let's calculate stress and strain at the 4 gauss points and then plot
% the deformed plate
Ue = U;
strainGP = zeros(3,4); % every point has 3 strain component e11, e22, e12
stressGP = zeros(3,4);  % every point has 3 stress component sigma_11. Sigma_22, Tau_12
gp = 0;

for g1 = 1:2
    for g2 = 1:2
        gp = gp + 1;
        xi = gaussian_points(g1);
        eta = gaussian_points(g2);

        weight = gaussian_weights(g1)*gaussian_weights(g2);
        % N = 0.25[(1-xi)(1-eta), (1+xi)(1-eta), (1+xi)(1+eta), (1-xi)(1+eta)
        % Row 1: Derivatives w.r.t xi
        % Row 2: Derivatives w.r.t eta
        dN_natural = 0.25* [-(1-eta), (1-eta), (1+eta), -(1+eta);
            -(1-xi),  -(1+xi),  (1+xi), (1-xi)];

        % Jacobian Convention:
        % J = [dx/dxi   dy/dxi;
        %      dx/deta  dy/deta];

        J = dN_natural * xe;
        detJ = det(J);

        if detJ <= 0
            error("Jacobian determinant is non-positive");
        end

        % Row1: dN/dx
        % Row2: dN/dy

        dN_physical = J \ dN_natural;

        % Strain displacement Matrix

        B = zeros(3, 8);

        for a = 1:4
            dNdx = dN_physical(1,a);
            dNdy = dN_physical(2,a);
            B(:, 2*a-1:2*a) = [dNdx 0;
                0 dNdy;
                dNdy dNdx];
        end
        strainGP(:, gp) = B*Ue;
        stressGP(:, gp) = D*strainGP(:, gp);
    end
end

disp('Each column corresponds to one Gauss point');
disp('Strains: rows = [epsilon_xx, epsilon_yy, gamma_xy]');
disp(strainGP);
disp('Stress: rows = [sigma_xx, sigma_yy, tau_xy]');
disp(stressGP/1e6);

% Now it's time to plot the deformed shape
% Plot original and exaggerated deformed geometry

scale = 40000;
deformed = coords + scale*[u,v];
% Repeat node 1 to close the outline
outline = [1,2,3,4,1];

figure;
plot(coords(outline,1), coords(outline,2), LineStyle="-", Color="g", LineWidth=1.5, DisplayName='Original');
hold on;
plot(deformed(outline,1), deformed(outline,2), LineStyle="--", Color="r", LineWidth=1.5, DisplayName='Deformed');
axis equal;
grid on;
xlabel('x (mm)');
ylabel('y (mm)');
title(sprintf('Plate deformation- magnification = %g', scale));
legend(Location="best");
hold off;