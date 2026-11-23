close all;
clear;
clc;

%% 1. Geometry
Rhole = 0.1;       % Hole radius [m]
W = 0.5;          % Quarter-model width [m]
H = 0.5;          % Quarter-model height [m]

%% 2. Material and loading
E = 200e9;        % Young's modulus [Pa]
nu = 0.3;
t = 0.01;         % Thickness [m]
T = 1e6;          % Right-edge horizontal traction [Pa]

D = E/(1-nu^2) * ...
    [1,  nu, 0;
     nu, 1,  0;
     0,  0,  (1-nu)/2];

%% 3. Mesh settings
nr = 32;          % Radial elements
nt = 64;         % Angular elements

% Square geometry: even nt places a ray through the outer corner
assert(W == H && mod(nt,2) == 0, ...
    'This mesh setup requires W = H and even nt.');

assert(Rhole > 0 && Rhole < min(W,H), ...
    'Hole radius must be positive and smaller than W and H.');

thetaNodes = linspace(0, pi/2, nt+1);

nNode = (nr+1)*(nt+1);
nElem = nr*nt;

coords = zeros(nNode, 2);
elements = zeros(nElem, 4);

%% 4. Generate coordinates
for j = 1:nt+1
    theta = thetaNodes(j);

    c = cos(theta);
    s = sin(theta);

    rRight = inf;
    rTop = inf;

    if c > 1e-12
        rRight = W/c;
    end

    if s > 1e-12
        rTop = H/s;
    end

    rOuter = min(rRight, rTop);
    radii = linspace(Rhole, rOuter, nr+1);

    for i = 1:nr+1
        node = (j-1)*(nr+1) + i;
        coords(node,:) = [radii(i)*c, radii(i)*s];
    end
end

%% 5. Generate counterclockwise connectivity
for j = 1:nt
    for i = 1:nr
        e = (j-1)*nr + i;

        n1 = (j-1)*(nr+1) + i;
        n2 = n1 + 1;
        n4 = n1 + nr + 1;
        n3 = n4 + 1;

        elements(e,:) = [n1, n2, n3, n4];
    end
end

%% 6. Store global DOFs for every element
% Each row: [u1 v1 u2 v2 u3 v3 u4 v4]
elementDofs = zeros(nElem, 8);

for a = 1:4
    elementDofs(:,2*a-1) = 2*elements(:,a)-1;
    elementDofs(:,2*a)   = 2*elements(:,a);
end

%% 7. Identify boundary nodes
tol = 1e-10*max(W,H);

leftNodes   = find(abs(coords(:,1)) < tol);
bottomNodes = find(abs(coords(:,2)) < tol);
rightNodes  = find(abs(coords(:,1)-W) < tol);
topNodes    = find(abs(coords(:,2)-H) < tol);

holeNodes = (1:(nr+1):nNode).';

% Sort outer-boundary nodes for loading and plotting
[~, order] = sort(coords(rightNodes,2));
rightNodes = rightNodes(order);

[~, order] = sort(coords(topNodes,1));
topNodes = topNodes(order);

%% 8. Initialize global system
nDof = 2*nNode;

K = sparse(nDof, nDof);
F = zeros(nDof, 1);

gaussPoints = [-1, 1]/sqrt(3);
gaussWeights = [1, 1];

minDetJ = inf;

%% 9. Stiffness integration and assembly
for e = 1:nElem
    nodes = elements(e,:);
    xe = coords(nodes,:);
    dofs = elementDofs(e,:);

    ke = zeros(8,8);

    for g1 = 1:2
        for g2 = 1:2
            xi = gaussPoints(g1);
            eta = gaussPoints(g2);

            weight = gaussWeights(g1)*gaussWeights(g2);

            [~, B, ~, detJ] = q4Shape(xe, xi, eta);

            minDetJ = min(minDetJ, detJ);

            ke = ke + B.'*D*B*t*detJ*weight;
        end
    end

    K(dofs,dofs) = K(dofs,dofs) + ke;
end

%% 10. Right-edge traction
for k = 1:length(rightNodes)-1
    n1 = rightNodes(k);
    n2 = rightNodes(k+1);

    edgeLength = norm(coords(n2,:) - coords(n1,:));
    edgeForce = T*t*edgeLength;

    F(2*n1-1) = F(2*n1-1) + edgeForce/2;
    F(2*n2-1) = F(2*n2-1) + edgeForce/2;
end

% No applied loads on the hole or top boundary

%% 11. Symmetry conditions and solution
% Left symmetry boundary: u = 0
% Bottom symmetry boundary: v = 0
fixedDofs = unique([2*leftNodes-1;
                   2*bottomNodes]);

freeDofs = setdiff((1:nDof).', fixedDofs);

U = zeros(nDof,1);
U(freeDofs) = K(freeDofs,freeDofs) \ F(freeDofs);

u = U(1:2:end);
v = U(2:2:end);

reactions = K*U - F;

%% 12. Recover Gauss-point strains and stresses
% Rows: [epsilon_xx; epsilon_yy; gamma_xy]
strainGP = zeros(3,4,nElem);

% Rows: [sigma_xx; sigma_yy; tau_xy]
stressGP = zeros(3,4,nElem);

gpCoords = zeros(4,2,nElem);

for e = 1:nElem
    nodes = elements(e,:);
    xe = coords(nodes,:);
    Ue = U(elementDofs(e,:));

    gp = 0;

    for g1 = 1:2
        for g2 = 1:2
            gp = gp + 1;

            xi = gaussPoints(g1);
            eta = gaussPoints(g2);

            [N, B] = q4Shape(xe, xi, eta);

            strainGP(:,gp,e) = B*Ue;
            stressGP(:,gp,e) = D*strainGP(:,gp,e);
            gpCoords(gp,:,e) = N*xe;
        end
    end
end

%% 13. Element-average horizontal stress
sigmaXX = zeros(nElem,1);

for e = 1:nElem
    sigmaXX(e) = mean(stressGP(1,:,e));
end

%% 14. Locate the largest sampled horizontal stress
sigmaSamples = reshape(stressGP(1,:,:), 4, nElem);

[maxSigma, index] = max(sigmaSamples(:));
[gpMax, elemMax] = ind2sub(size(sigmaSamples), index);

peakPosition = gpCoords(gpMax,:,elemMax);

%% 15. Hole-boundary hoop stress
% Evaluate at the midpoint of each hole-facing straight edge
thetaBoundary = zeros(nt,1);
hoopStress = zeros(nt,1);

for j = 1:nt
    e = (j-1)*nr + 1;

    nodes = elements(e,:);
    xe = coords(nodes,:);
    Ue = U(elementDofs(e,:));

    % Hole-facing edge: xi = -1
    % Edge midpoint: eta = 0
    [N, B] = q4Shape(xe, -1, 0);

    stressBoundary = D*(B*Ue);
    position = N*xe;

    theta = atan2(position(2), position(1));
    thetaBoundary(j) = theta;

    sx = stressBoundary(1);
    sy = stressBoundary(2);
    txy = stressBoundary(3);

    hoopStress(j) = sx*sin(theta)^2 ...
                  + sy*cos(theta)^2 ...
                  - 2*txy*sin(theta)*cos(theta);
end

[maxHoop, peakIndex] = max(hoopStress);

%% 16. Print checks and results
fprintf('Nodes    = %d\n', nNode);
fprintf('Elements = %d\n', nElem);
fprintf('DOFs     = %d\n\n', nDof);

fprintf('Minimum Gauss-point det(J) = %.6e m^2\n', minDetJ);

fprintf('Applied horizontal force = %.6f N\n', ...
        sum(F(1:2:end)));

fprintf('Horizontal reaction = %.6f N\n', ...
        sum(reactions(1:2:end)));

fprintf('Vertical reaction = %.6f N\n', ...
        sum(reactions(2:2:end)));

fprintf('Maximum horizontal displacement = %.6f mm\n\n', ...
        max(u)*1000);

fprintf('Maximum Gauss-point sigma_xx = %.6f MPa\n', ...
        maxSigma/1e6);

fprintf('Minimum Gauss-point sigma_xx = %.6f MPa\n', ...
        min(sigmaSamples(:))/1e6);

fprintf('Maximum sampled sigma_xx/T = %.6f\n', ...
        maxSigma/T);

fprintf('Peak position: x = %.6f m, y = %.6f m\n\n', ...
        peakPosition(1), peakPosition(2));

fprintf('Maximum sampled boundary hoop stress = %.6f MPa\n', ...
        maxHoop/1e6);

fprintf('Maximum sampled boundary hoop stress/T = %.6f\n', ...
        maxHoop/T);

fprintf('Boundary peak angle = %.6f degrees\n', ...
        thetaBoundary(peakIndex)*180/pi);

%% 17. Plot mesh and boundaries
figure;

patch('Faces', elements, 'Vertices', coords, ...
      'FaceColor', 'none', 'EdgeColor', [0.6,0.6,0.6], ...
      'HandleVisibility', 'off');

hold on;

plot(coords(holeNodes,1), coords(holeNodes,2), 'r-', ...
    'LineWidth', 2, 'DisplayName', 'Hole');

plot(coords(leftNodes,1), coords(leftNodes,2), 'b-', ...
    'LineWidth', 2, 'DisplayName', 'Left symmetry');

plot(coords(bottomNodes,1), coords(bottomNodes,2), 'g-', ...
    'LineWidth', 2, 'DisplayName', 'Bottom symmetry');

plot(coords(rightNodes,1), coords(rightNodes,2), 'm-', ...
    'LineWidth', 2, 'DisplayName', 'Loaded edge');

plot(coords(topNodes,1), coords(topNodes,2), 'c-', ...
    'LineWidth', 2, 'DisplayName', 'Free top edge');

axis equal;
xlabel('x [m]');
ylabel('y [m]');
title('Quarter-plate Q4 mesh');
legend('Location', 'bestoutside');
grid on;
hold off;

%% 18. Plot deformation
scale = 1000;

originalMM = coords*1000;
deformedMM = (coords + scale*[u,v])*1000;

figure;

hOriginal = patch('Faces', elements, 'Vertices', originalMM, ...
    'FaceColor', 'none', 'EdgeColor', [0.6,0.6,0.6]);

hold on;

hDeformed = patch('Faces', elements, 'Vertices', deformedMM, ...
    'FaceColor', 'b', 'FaceAlpha', 0.15, ...
    'EdgeColor', 'b');

axis equal;
xlabel('x [mm]');
ylabel('y [mm]');
title(sprintf('Deformation magnification = %g', scale));
legend([hOriginal,hDeformed], {'Original','Deformed'});
grid on;
hold off;

%% 19. Plot element-average horizontal stress
figure;

patch('Faces', elements, 'Vertices', originalMM, ...
      'FaceVertexCData', sigmaXX/1e6, ...
      'FaceColor', 'flat', 'EdgeColor', 'none');

axis equal;
xlabel('x [mm]');
ylabel('y [mm]');
title('Element-average sigma_{xx} [MPa]');
colorbar;
grid on;

%% 20. Plot hoop stress against angle
thetaReference = linspace(0, pi/2, 200);

% Infinite-plate reference, not an exact finite-plate solution
hoopReference = T*(1 - 2*cos(2*thetaReference));

figure;

plot(thetaBoundary*180/pi, hoopStress/T, 'b-', ...
    'LineWidth', 1.5, 'DisplayName', 'Finite-plate FEM');

hold on;

plot(thetaReference*180/pi, hoopReference/T, 'r--', ...
    'LineWidth', 1.5, 'DisplayName', 'Infinite-plate reference');

xlabel('Angle from positive x-axis [degrees]');
ylabel('Hoop stress / applied traction');
title('Hoop stress around the hole');
legend('Location', 'best');
grid on;
hold off;