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
Plate_coords = [0 0;  % node 1 : bottom left
                L 0;  % node 2: bottom right
                L H;  % node 3: top right
                0 H]; % node 4: top left
% structured mesh
nx = 4; % Element along x
ny = 4; % Element along y

xNodes = linspace(0,L,nx+1);
yNodes = linspace(0,H,ny+1);

nNode = (nx+1)*(ny+1);
coords = zeros(nNode, 2);

% Number nodes from left to right, then bottom to top

for j= 1:ny+1
    for i = 1:nx+1
        node = (j-1)*(nx+1)+i;
        coords(node,:) = [xNodes(i) yNodes(j)];
    end
end

% counter-clockwise connectivity
% Local order: bottom-left, bottom-right, top-right, top-left

nElem = nx*ny;
elements = zeros(nElem, 4);

for j=1:ny
    for i= 1:nx
        e=(j-1)*nx+i;
        n1 = (j-1)*(nx+1) +i;
        n2 = n1 + 1;
        n4 = n1 + nx + 1;
        n3 = n4 + 1;
        elements(e,:) = [n1 n2 n3 n4];
    end
end

% Degree of freedom
% ordering: [u1, v1, u2, v2, u3, v3, u4, v4,....]

nDof = 2 * nNode;

K = sparse(nDof, nDof);
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

% Element stiffness integration and assembly

for e= 1:nElem
    nodes = elements(e,:);
    xe = coords(nodes, :);
    
    % convert 4 node numbers into 8 DOF numbers

    dofs = zeros(1,8);

    for a = 1:4
        dofs(2*a-1 : 2*a) = [2*nodes(a)-1. 2*nodes(a)];
    end

    ke = zeros(8,8);

    for g1 = 1:2
        for g2 = 1:2
            xi = gaussian_points(g1);
            eta = gaussian_points(g2);
            weight = gaussian_weights(g1) * gaussian_weights(g2);

            % Natural Derivatives: rows corresponds to xi and eta

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
                error('Jacobian determinant is non-positivein element %d',e);
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

    K(dofs, dofs) = K(dofs, dofs) + ke;
end

% Assemble uniform traction on the right boundary
% Each straight boundary segment contributes half its force to each of its two endpoints

for j = 1:ny
    bottomNode = (j-1)*(nx+1) + nx +1;
    topNode = bottomNode + nx + 1;

    edgelength = norm(coords(topNode, :) - coords(bottomNode, :));
    edgeForce = T*t*edgelength;

    F(2*bottomNode-1) = F(2*bottomNode-1) + edgeForce/2;
    F(2*topNode-1) = F(2*topNode-1) + edgeForce/2;
end

% Apply displacement boundary conditions
% left edge : horizontal displacement = 0;
%Bottom-left corner: vertical displacement = 0
% other vertical displacement remain free

tol = 1e-10 * max(L,H);

leftNodes = find(abs(coords(:,1))<tol);
bottomLeftNodes = 1;

fixedDofs = unique([2*leftNodes-1 ; 2* bottomLeftNodes]);

alldofs = (1:nDof).';
freedofs = setdiff(alldofs, fixedDofs);


% Solve 
U = zeros(nDof,1);
U(freedofs) = K(freedofs, freedofs) \ F(freedofs);

% Extract displacement components:

u = U(1:2:end);
v = U(2:2:end);

R = K*U -F;

% Analytical Displacement Comparison
% Exact solution for the chosen loading and constraints

uExact = T/E *coords(:,1);
vExact = -nu*T/E *coords(:,2);

fprintf('Nomber of nodes        = %d\n',nNode);
fprintf('Number of elements     = %d\n',nElem);
fprintf('Number of DOFs         = %d\n',nDof);
fprintf('Number of fixed DOFs   = %d\n',numel(fixedDofs));
fprintf('Number of free DOFs    = %d\n',numel(freedofs));

%disp('Node          u,mm        v,mm');
%disp([(1:nNode).', u*1000, v*1000]);

fprintf('Total horizontal force applied = %.6f N\n', sum(F(1:2:end)));
fprintf('Total vertical force applied = %.6f N\n', sum(F(2:2:end)));

fprintf('Total Horizontal Reaction = %.6f N\n', sum(R(1:2:end)));
fprintf('Total Vertical Reaction = %.6f N\n', sum(R(2:2:end)));

% Display nodal displacement errors

%disp('Nodal Displacement Error');
%disp('  Node       u_Error(%)       v_Error(%)');
%disp([(1:nNode).', abs((u-uExact)*1000), abs((v-vExact)*1000)]);


% Now let's calculate stress and strain at all gauss points and then plot the deformed plate
% Rows: Strain/Stress components
% Columns: Gauss points
% Third index: element number

strainGP = zeros(3,4,nElem); % every point has 3 strain component e11, e22, e12
stressGP = zeros(3,4,nElem);  % every point has 3 stress component sigma_11. Sigma_22, Tau_12

for e = 1:nElem
    nodes = elements(e, :);
    xe = coords(nodes, :);

    dofs = zeros(1,8);

    for a = 1:4
        dofs(2*a-1:2*a) = [2*nodes(a)-1, 2*nodes(a)];
    end

    Ue = U(dofs);
    gp =0;

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
            strainGP(:, gp, e) = B*Ue;
            stressGP(:, gp, e) = D*strainGP(:, gp,e);
        end
    end
end

% Mean Gauss-point sigma_xx for each element

%sigma_xx = squeeze(mean(stressGP(1,:,:), 2));
sigma_xx = zeros(nElem, 1);

for e = 1:nElem
    sigma_xx(e) = mean(stressGP(1,:,e));
end

fprintf('\nMinimum Gauss-point sigma_xx = %0.6f Mpa\n', min(reshape(stressGP(1,:,:), [], 1))/1e6);
fprintf('Maximum Gauss-point sigma_xx = %0.6f Mpa\n', max(reshape(stressGP(1,:,:), [], 1))/1e6);
fprintf('Mean Gauss-point sigma_xx = %0.6f Mpa\n', mean(sigma_xx)/1e6);


% Plot mesh and node numbers

figure;
patch(Vertices=coords, Faces=elements, FaceColor="cyan", EdgeColor="b", LineWidth=1);
axis equal;
grid on;
xlabel('x (mm)');
ylabel('y (mm)');
hold on;

for node = 1:nNode
    text(coords(node,1), coords(node,2), sprintf(' %d', node), Color="r", FontSize=8);
end
title('Q4 Mesh with global node numbers');
hold off;


% Now it's time to plot the deformed shape
% Plot original and exaggerated deformed geometry

scale = 50000;
original_mm = coords*1000;
deformed_mm = (coords + scale*[u,v])*1000;
% Repeat node 1 to close the outline

figure;
hOriginal = patch( "Faces",elements, "Vertices",original_mm, FaceColor = 'g', FaceAlpha= 0.5, EdgeColor = 'yellow');
hold on;
hDeformed = patch( "Faces",elements, "Vertices",deformed_mm, FaceColor = 'b', facealpha =0.2, EdgeColor = 'blue' );
axis equal;
title(sprintf('Plate deformation- magnification = %g', scale));
grid on;
xlabel('x (mm)');
ylabel('y (mm)');
legend([hOriginal , hDeformed], {'Original','Deformed'}, Location="best");
hold off;

% PLot element- average horizontal stress
figure;
patch(Vertices=original_mm, Faces=elements, FaceColor="flat", FaceVertexCData=sigma_xx/1e6, EdgeColor="none", FaceAlpha = 1);
axis equal;
grid on;
xlabel('x (mm)');
ylabel('y (mm)');
colorbar;
title('Element- average sigma_xx, MPa');

