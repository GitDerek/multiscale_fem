function Untitled()
tic

%% define a uniform finite element mesh over the unit square [0,1]x[0,1]
epsilon = 0.001;
a1 = 1/2;
a2 = sqrt(3)/2;

n = 20;
Th = RectangleMesh(0,0,1,1,n,n);
coordinates = Th.coordinates;
elements = Th.elements;

%% compute boundary vertices and edges
% four boundary edges
k1 = find( coordinates(:,2)==0  );
k2 = find( coordinates(:,1)==1 & coordinates(:,2)<=1   & coordinates(:,2) >=0);
k3 = find( coordinates(:,2)==1  );
k4 = find( coordinates(:,1)==0 & coordinates(:,2)<=1   & coordinates(:,2) >=0);
%
dirichlet = sort(unique([k1;k2;k3;k4]));
%
nE = size(elements,1);                                          % number of elemenet
nC = size(coordinates,1);                                            % number of total vertices

freenodes = setdiff(1:nC,dirichlet)';

%% perform assembly of global Laplacian matrix
tic
I = reshape(elements(:, [1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4])', 16 * nE, 1);
J = reshape(elements(:, [1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4])', 16 * nE, 1);
A = zeros(16*nE, 1);
for i = 1:nE
    node = Th.elements(i,:);
    idx = 16*(i-1)+1:16*i;
    A(idx) = Q1elem_Lap(Th.coordinates(node,:), epsilon, a1, a2);
end
A = sparse(I,J,A,nC,nC);
toc

%% perform assembly of global souce verctor
rhs = zeros(nC,1);

%%
uh = zeros(nC,1);
uh(coordinates(:,1)<=0.5 & coordinates(:, 2) == 0) = 1;
uh(k4) = 1;
uh(k3) = 0;
uh(k2) = 0;
rhs = rhs - A * uh;
uh(freenodes) = A(freenodes,freenodes)\rhs(freenodes);

toc
%%
close all
figure(1)
trimesh(elements, coordinates(:, 1), coordinates(:, 2), uh)
set(gca,'fontsize',16)
title (['approximate solution: h = 1/',num2str(n)]);
xlabel('x');
ylabel('y');
end

function M = Q1elem_Lap(vertices, epsilon, a1, a2)

node = vertices;
[Th, macroMatrix, msBasis] = solveLocalBasis(node(1,1), node(1,2), node(3,1), node(3,2), 20, epsilon, a1, a2);
phi = Q1basisCoefficient(node(1,1), node(3,1), node(1,2), node(3,2), Th.coordinates);
M = msBasis' * macroMatrix * msBasis;

end

function M = assumbleLocalMatrix(vertices, epsilon, a1, a2)
% computation of the element stiffness matrix with Q1 element
% for Laplace operator

node = vertices;
M = zeros(4,4);

% set up 2x2 Gauss points
gpt = 1.0/sqrt(3.0);
s(1) = -gpt;  t(1) = -gpt;
s(2) =  gpt;  t(2) = -gpt;
s(3) =  gpt;  t(3) =  gpt;
s(4) = -gpt;  t(4) =  gpt;


for igpt = 1:4               %  start gauss quadrature
    sigpt = s(igpt);
    tigpt = t(igpt);
    
    %  evaluate derivatives etc
    [phi, dphids, dphidt] = Q1basis(sigpt, tigpt);
    
    dxds = node(1,1) * dphids(1) + node(2,1) * dphids(2) + node(3,1) * dphids(3)+node(4,1)*dphids(4);
    dxdt = node(1,1) * dphidt(1) + node(2,1) * dphidt(2) + node(3,1) * dphidt(3)+node(4,1)*dphidt(4);
    dyds = node(1,2) * dphids(1) + node(2,2) * dphids(2) + node(3,2) * dphids(3)+node(4,2)*dphids(4);
    dydt = node(1,2) * dphidt(1) + node(2,2) * dphidt(2) + node(3,2) * dphidt(3)+node(4,2)*dphidt(4);
    
    % compute Jacobian
    jac = dxds*dydt-dxdt*dyds;
    
    
    % compute derivate of shape function
    for ivtx = 1:4
        dphidx(ivtx) = ( dphids(ivtx)*dydt - dphidt(ivtx)*dyds) ./ jac;
        dphidy(ivtx) = (-dphids(ivtx)*dxdt + dphidt(ivtx)*dxds) ./ jac;
    end
    
    % Local sitiffness matrix
    for j = 1:4
        for i = 1:4
            M(i,j) = M(i,j) + epsilon .* dphidx(i) * dphidx(j) * jac;
            M(i,j) = M(i,j) + epsilon .* dphidy(i) * dphidy(j) * jac;
            M(i,j) = M(i,j) + (a1 .* dphidx(j) .* phi(i) + a2 .* dphidy(j) .* phi(i)) *jac;
        end
    end
    
end % end of gauss quadrature
end

function [Th, A, msBasis] = solveLocalBasis(x0, y0, x1, y1, n, epsilon, a1, a2)

Th = RectangleMesh(x0, y0, x1, y1, n, n);
macroCoordinates = Th.coordinates;
macroElements = Th.elements;

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

freenodes = setdiff(1:nC,dirichlet)';

% %% perform assembly of global Laplacian matrix
I = reshape(macroElements(:, [1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4])', 16 * nE, 1);
J = reshape(macroElements(:, [1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4])', 16 * nE, 1);
A = zeros(16*nE, 1);
for i = 1:nE
    node = macroElements(i,:);
    idx = 16*(i-1)+1:16*i;
    A(idx) = assumbleLocalMatrix(macroCoordinates(node,:), epsilon, a1, a2);
end
A = sparse(I,J,A,nC,nC);

%% perform assembly of global souce verctor
msBasis = zeros(nC ,4);
for basis_no = 1:4
    msBasis(dirichlet,basis_no) = macroElementBoundaryCondition(x0, x1, y0, y1, macroCoordinates(dirichlet,1), macroCoordinates(dirichlet,2), basis_no);
    rhs = zeros(nC,1);
    rhs = rhs - A * msBasis(:, basis_no);
    msBasis(freenodes, basis_no) = A(freenodes,freenodes)\rhs(freenodes);
end

a = 1 - (msBasis(:, 1) + msBasis(:, 3) + msBasis(:, 3));

close all;
figure(99)
% for basis_no = 1:4
subplot(2,2,basis_no)
trimesh(macroElements, macroCoordinates(:, 1), macroCoordinates(:, 2), a)
% set(gca,'fontsize',16)
title (['Multiscale Basis ', num2str(basis_no)]);
xlabel('x');
ylabel('y');
axis([x0, x1, y0, y1, 0,1])
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

function phi = Q1basisCoefficient(x0, x1, y0, y1, coordinates)

hx = x1 - x0;
hy = y1 - y0;
s = coordinates(:,1);
t = coordinates(:,2);
nc = size(coordinates,1);

phi = zeros(nc, 4);

phi(:, 1) = (1 - (s - x0) ./hx) .* (1 - (t - y0) ./ hy);
phi(:, 2) = (s - x0) ./hx .* (1 - (t - y0) ./ hy);
phi(:, 3) = (s - x0) ./hx .* (t - y0) ./ hy;
phi(:, 4) = (1 - (s - x0) ./hx) .* (t - y0) ./ hy;

end

function [phi,dphids,dphidt] = Q1basis(s,t)

% definition of Q1 shape function
phi(1) = 0.25*(s-1.)*(t-1.);
phi(2) = -0.25*(s+1.)*(t-1.);
phi(3) = 0.25*(s+1.)*(t+1.);
phi(4) = -0.25*(s-1.)*(t+1.);
%
dphids(1) = 0.25*(t-1.);
dphids(2) = -0.25*(t-1.);
dphids(3) = 0.25*(t+1.);
dphids(4) = -0.25*(t+1.);
%
dphidt(1) = 0.25*(s-1.);
dphidt(2) = -0.25*(s+1.);
dphidt(3) = 0.25*(s+1.);
dphidt(4) = -0.25*(s-1.);
end