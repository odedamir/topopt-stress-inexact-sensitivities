%% FUNCTION mgcg_adj - MULTIGRID PRECONDITIONED CONJUGATE GRADIENTS
function [i,relres,u] = mgcg_adj(A,b,u,Lfac,Ufac,Pu,nl,nswp,tol,maxiter,...
    u_state,edofMat,KE,x,penal,Emax,Emin)
r = b - A{1,1}*u;
res0 = norm(b);
% Alternative convergence
lamdku1(1:size(x),1) = -1000;
lamdku2(1:size(x),1) = 1000;
% Jacobi smoother
omega = 0.8;
invD = cell(nl-1,1);
for l = 1:nl-1
    invD{l,1}= 1./spdiags(A{l,1},0);
end
for i = 1:1e6 
    z = VCycle(A,r,Lfac,Ufac,Pu,1,nl,invD,omega,nswp);
    rho = r'*z;
    if i == 1
        p = z;
    else
        beta = rho/rho_p;
        p = beta*p + z;
    end
    q = A{1,1}*p;
    dpr = p'*q;
    alpha = rho/dpr;
    u = u + alpha*p;
    r = r - alpha*q;
    rho_p = rho;
    relres = norm(r)/res0;
    ce1 = sum(u(edofMat)*KE.*u_state(edofMat),2); % Lam'*K*U
    lamdku = penal*(Emax-Emin)*x.^(penal-1).*ce1; % Lam'*dKdrho*U
    altres1 = (lamdku-lamdku2)/2;
    altres2 = (lamdku-lamdku1);
    scale = max(abs(lamdku));
    altres = max(max(abs(altres1),abs(altres2)))/scale;
%     fprintf(' MGCGiter:%3i relres:%4.3f altres:%4.3f \n',...
%         i,relres,altres);
    if altres < tol || i>=maxiter
        break
    end
    lamdku2 = lamdku1;
    lamdku1 = lamdku;
end
end