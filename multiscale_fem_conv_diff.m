clear all, format long
%% define a uniform finite element mesh over the unit square [0,1]x[0,1]
epsilon = 0.001;
a1 = 1/2;
a2 = sqrt(3)/2;

n = 8;
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
hE = sqrt(2)*(1/n);                                         % ;
hF = (1/n);
freenodes = setdiff(1:nC,dirichlet)';

%% perform assembly of global Laplacian matrix
diffusionGlobalMatrix = sparse(nC,nC);
for k = 1:nE
    diffusionGlobalMatrix(elements(k,:),elements(k,:)) ...
        = diffusionGlobalMatrix(elements(k,:),elements(k,:)) + q1elem_Lap(coordinates(elements(k,:),:), epsilon ,a1, a2);
end

%% perform assembly of global souce verctor
rhs = zeros(nC,1);

%%
uh = zeros(nC,1);
uh(coordinates(:,1)<=0.5 & coordinates(:, 2) == 0) = 1;
uh(k4) = 1;
uh(k3) = 0;
uh(k2) = 0;
rhs = rhs - diffusionGlobalMatrix * uh;
uh(freenodes) = diffusionGlobalMatrix(freenodes,freenodes)\rhs(freenodes);

%%
close all
figure(1)
trimesh(elements, coordinates(:, 1), coordinates(:, 2), uh)
set(gca,'fontsize',16)
title (['approximate solution: h = 1/',num2str(n)]);
xlabel('x');
ylabel('y');