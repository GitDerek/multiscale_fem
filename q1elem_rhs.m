function elementsRHS = Q1elem_rhs(vertices)

node = vertices;
elementsRHS = zeros(4,1);

% set up 2x2 Gauss points
gpt=1.0/sqrt(3.0);
s(1) = -gpt;  t(1) = -gpt;
s(2) =  gpt;  t(2) = -gpt;
s(3) =  gpt;  t(3) =  gpt;
s(4) = -gpt;  t(4) =  gpt;

for igpt = 1:4               %  start gauss quadrature
    sigpt=s(igpt);
    tigpt=t(igpt);
    
    %  evaluate derivatives etc
    [phi,dphids,dphidt] = Q1basis(sigpt,tigpt);
    
    dxds = node(1,1)*dphids(1)+node(2,1)*dphids(2)+node(3,1)*dphids(3)+node(4,1)*dphids(4);
    dxdt = node(1,1)*dphidt(1)+node(2,1)*dphidt(2)+node(3,1)*dphidt(3)+node(4,1)*dphidt(4);
    dyds = node(1,2)*dphids(1)+node(2,2)*dphids(2)+node(3,2)*dphids(3)+node(4,2)*dphids(4);
    dydt = node(1,2)*dphidt(1)+node(2,2)*dphidt(2)+node(3,2)*dphidt(3)+node(4,2)*dphidt(4);
    
    % jacobian
    jac = dxds*dydt-dxdt*dyds;
    
    for ivtx = 1:4
        dphidx(ivtx) = ( dphids(ivtx)*dydt - dphidt(ivtx)*dyds)./jac;
        dphidy(ivtx) = (-dphids(ivtx)*dxdt + dphidt(ivtx)*dxds)./jac;
    end
    
    %
    xx = node(1,1)*phi(1)+node(2,1)*phi(2)+node(3,1)*phi(3)+node(4,1)*phi(4);
    yy = node(1,2)*phi(1)+node(2,2)*phi(2)+node(3,2)*phi(3)+node(4,2)*phi(4);
    
    elementsSourceTerm = SourceTerm(xx,yy);
    for j = 1:4
        elementsRHS(j) = elementsRHS(j) + elementsSourceTerm*phi(j)*jac;
    end
end

