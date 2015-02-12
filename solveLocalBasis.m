function [Th, msBasis] = solveLocalBasis(x0, y0, x1, y1, n, epsilon, a1, a2)

Th = RectangleMesh(x0, y0, x1, y1, n, n);
macroCoordinates = Th.coordinates;
macroElements = Th.elements;

% Th.showMesh()

%% compute boundary vertices and edges
% four boundary edges
k1 = find( macroCoordinates(:,2) == y0 );
k2 = find( macroCoordinates(:,1) == x1 );
k3 = find( macroCoordinates(:,2) == y1 );
k4 = find( macroCoordinates(:,1) == x0 );
%
dirichlet = sort(unique([k1;k2;k3;k4]));
%
nE = Th.getNumElements();                                          % number of elemenet
nC = Th.getNumCoordinates();                                            % number of total vertices
hE = sqrt(2)*(1/n);                                         % ;
hF = (1/n);
freenodes = setdiff(1:nC,dirichlet)';

% %% perform assembly of global Laplacian matrix
localBasis = sparse(nC,nC);
for k = 1:nE
    localBasis(macroElements(k,:),macroElements(k,:)) ...
        = localBasis(macroElements(k,:),macroElements(k,:)) + assumbleLocalMatrix(macroCoordinates(macroElements(k,:),:), epsilon, a1, a2);
end

%% perform assembly of global souce verctor
msBasis = zeros(nC ,4);
for basis_no = 1:4
    msBasis(dirichlet,basis_no) = macroElementBoundaryCondition(x0, x1, y0, y1, macroCoordinates(dirichlet,1), macroCoordinates(dirichlet,2), basis_no);
    rhs = zeros(nC,1);
    rhs = rhs - localBasis * msBasis(:, basis_no);
    msBasis(freenodes, basis_no) = localBasis(freenodes,freenodes)\rhs(freenodes);
end

% close all;
% figure(1)
% for basis_no = 1:4
% subplot(2,2,basis_no)
% trimesh(macroElements, macroCoordinates(:, 1), macroCoordinates(:, 2), msBasis(:, basis_no),...
%     'facecolor', 'interp', 'linestyle', 'none')
% set(gca,'fontsize',16)
% title (['Multiscale Basis ', num2str(basis_no)]);
% xlabel('x');
% ylabel('y');
% axis([x0, x1, y0, y1, 0,1])
% end
end

function phi = macroElementBoundaryCondition(x0, x1, y0, y1, s, t, basis_no)

hx = x1 - x0;
hy = y1 - y0;
switch basis_no
    case 1
        phi = (1 - (s - x0) ./hx) .* (1 - (t - y0) ./ hy);
    case 2
        phi = (s - x0) ./hx .* (1 - (t - y0) ./ hy);
    case 3
        phi = (s - x0) ./hx .* (t - y0) ./ hy;
    case 4
        phi = (1 - (s - x0) ./hx) .* (t - y0) ./ hy;
    otherwise
        disp('Unknown boundary condition')
end
end