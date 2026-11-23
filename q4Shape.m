function [N, B, J, detJ] = q4Shape(xe, xi, eta)
% xe: 4-by-2 element coordinates
% xi, eta: natural coordinates of the evaluation point

N = 0.25 * ...
    [(1-xi)*(1-eta), ...
    (1+xi)*(1-eta), ...
    (1+xi)*(1+eta), ...
    (1-xi)*(1+eta)];

dNnatural = 0.25 * ...
    [-(1-eta),  (1-eta), (1+eta), -(1+eta);
    -(1-xi), -(1+xi),  (1+xi),   (1-xi)];

J = dNnatural*xe;
detJ = det(J);

if detJ <= 0
    error('Nonpositive Jacobian determinant.');
end

dNphysical = J \ dNnatural;

B = zeros(3, 8);

for a = 1:4
    dNdx = dNphysical(1, a);
    dNdy = dNphysical(2, a);

    B(:, 2*a-1:2*a) = ...
        [dNdx, 0;
        0,    dNdy;
        dNdy, dNdx];
end
end