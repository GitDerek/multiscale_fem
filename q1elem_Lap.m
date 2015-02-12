function M = Q1elem_Lap(vertices, epsilon, a1, a2)
% computation of the element stiffness matrix with Q1 element
% for Laplace operator

node = vertices;
[Th, msBasis] = solveLocalBasis(node(1,1), node(1,2), node(3,1), node(3,2), 20, epsilon, a1, a2);

nE = Th.getNumElements();                                          % number of elemenet
nC = Th.getNumCoordinates();                                       

macroMatrix = sparse(nC,nC);
for k = 1:nE
    macroMatrix(Th.elements(k,:), Th.elements(k,:)) ...
        = macroMatrix(Th.elements(k,:), Th.elements(k,:)) + assumblelocal(Th.coordinates(Th.elements(k,:),:), epsilon, a1, a2);
end


M = msBasis' * macroMatrix * msBasis;

end

function M = assumblelocal(vertices, epsilon, a1, a2)
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
    [phi,dphids,dphidt] = q1basis(sigpt,tigpt);
    
    dxds = node(1,1)*dphids(1)+node(2,1)*dphids(2)+node(3,1)*dphids(3)+node(4,1)*dphids(4);
    dxdt = node(1,1)*dphidt(1)+node(2,1)*dphidt(2)+node(3,1)*dphidt(3)+node(4,1)*dphidt(4);
    dyds = node(1,2)*dphids(1)+node(2,2)*dphids(2)+node(3,2)*dphids(3)+node(4,2)*dphids(4);
    dydt = node(1,2)*dphidt(1)+node(2,2)*dphidt(2)+node(3,2)*dphidt(3)+node(4,2)*dphidt(4);
    
    % jacobian
    jac = dxds*dydt-dxdt*dyds;
    
    
    for ivtx = 1:4
        dphidx(ivtx) = ( dphids(ivtx)*dydt- dphidt(ivtx)*dyds)./jac;
        dphidy(ivtx) = (-dphids(ivtx)*dxdt+ dphidt(ivtx)*dxds)./jac;
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