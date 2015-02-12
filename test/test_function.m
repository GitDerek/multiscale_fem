function test_function

n = 4;
Th = RectangleMesh(0,0,1,1,n,n);
coordinates = Th.coordinates;
elements = Th.elements;

figure(2)
Th.showMesh();

cc = coordinates(elements(6,:),:)
solveLocalBasis(cc(1,1), cc(1,2), cc(3,1), cc(3,2), 100, 0.01, 1/2, sqrt(3)/2);
% [bdNode,bdEdge,isBdNode,isBdElem] = findquadboundary(Th.elements);
% [edge, elem2edge] = auxstructure(Th.elements)
end

function [bdNode,bdEdge,isBdNode,isBdElem] = findquadboundary(elem)
%% FINDQUADBOUNDARY finds the boundary of a quad mesh
%
%   [bdNode,bdEdge,isBdNode] = findquadboundary(elem) finds boundary nodes and
%   edges for a 2-dimensional quad mesh. Note only the topological structure of
%   the mesh is needed.
%
%   See also 
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

[NT, NV] = size(elem);
totalEdge = zeros(NV*NT,2);
totalEdge(:,1) = elem(:);
totalEdge(:,2) = elem([NT+1:NV*NT,1:NT]');
totalEdge(:) = sort(totalEdge,2);
[i,j,s] = find(sparse(double(totalEdge(:,1)),double(totalEdge(:,2)),1));
bdEdge = [i(s==1),j(s==1)];
isBdNode = false(max(elem(:)),1); 
isBdNode(bdEdge) = true;
bdNode = find(isBdNode);
isBdElem = isBdNode(elem(:,1)) | isBdNode(elem(:,2)) | isBdNode(elem(:,3)) | isBdNode(elem(:,4));
end

function [edge, elem2edge] = auxstructure(elem)
%% AUXSTRUCTURE auxiliary structure for a 2-D triangulation.
%
%  T = AUXSTRUCTURE(elem) constucts the indices map between elements, edges 
%  and nodes, and the boundary information. T is a structure. 
%
%  T.neighbor(1:NT,1:3): the indices map of neighbor information of elements, 
%  where neighbor(t,i) is the global index of the element oppoiste to the 
%  i-th vertex of the t-th element. 
%
%  T.elem2edge(1:NT,1:3): the indices map from elements to edges, elem2edge(t,i) 
%  is the edge opposite to the i-th vertex of the t-th element.
%
%  T.edge(1:NE,1:2): all edges, where edge(e,i) is the global index of the 
%  i-th vertex of the e-th edge, and edge(e,1) < edge(e,2) 
%
%  T.bdEdge(1:Nbd,1:2): boundary edges with positive oritentation, where
%  bdEdge(e,i) is the global index of the i-th vertex of the e-th edge for
%  i=1,2. The positive oritentation means that the interior of the domain
%  is on the left moving from bdEdge(e,1) to bdEdge(e,2). Note that this
%  requires elem is positive ordered, i.e., the signed area of each
%  triangle is positive. If not, use elem = fixorder(node,elem) to fix the
%  order.
%
%  T.edge2elem(1:NE,1:4): the indices map from edge to element, where 
%  edge2elem(e,1:2) are the global indexes of two elements sharing the e-th
%  edge, and edge2elem(e,3:4) are the local indices of e to edge2elem(e,1:2).
%
%  To save space all the data type in T is uint32. When use them as a input
%  of sparse(i,j,s,m,n), please change them into double type.
% 
%  See also auxstructure3.
% 
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

totalEdge = uint32(sort([elem(:,[4,1]); ...
                         elem(:,[1,2]); ...
                         elem(:,[2,3]); ...
                         elem(:,[3,4])], 2));

matlabversion = version;
if str2double(matlabversion(end-5:end-2)) > 2012
    [edge, i2, j] = unique(totalEdge,'rows','legacy');
else
    [edge, i2, j] = unique(totalEdge,'rows');
end
NT = size(elem,1);
elem2edge = uint32(reshape(j,NT,4));

% i1(j(3*NT:-1:1)) = 3*NT:-1:1; 
% i1 = i1';
% k1 = ceil(i1/NT); 
% k2 = ceil(i2/NT); 
% t1 = i1 - NT*(k1-1);
% t2 = i2 - NT*(k2-1);
% ix = (i1 ~= i2); 
% neighbor = uint32(accumarray([[t1(ix),k1(ix)];[t2,k2]],[t2(ix);t1],[NT 3]));
% edge2elem = uint32([t1,t2,k1,k2]);
% bdElem = t1(t1 == t2);
% bdk1 = k1(t1 == t2);
% bdEdge = [elem(bdElem(bdk1==1),[2 3]); elem(bdElem(bdk1==2),[3 1]);...
%           elem(bdElem(bdk1==3),[1 2])];
% T = struct('neighbor',neighbor,'elem2edge',elem2edge,'edge',edge,...
%            'edge2elem',edge2elem,'bdElem',bdElem,'bdEdge',bdEdge);
       
end