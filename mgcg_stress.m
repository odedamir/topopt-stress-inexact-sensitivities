%% FUNCTION mgcg_stress - MULTIGRID PRECONDITIONED CONJUGATE GRADIENTS WITH STRESS COMPUTATION
% If actual max and location are needed, use the commented function and
% uncomment some parts of the code
% function [i,relres,u,where,maxVM] = mgcg_stress(A,b,u,Lfac,Ufac,Pu,nl,nswp,tol,maxiter,stress_data)
function [i,relres,u] = mgcg_stress(A,b,u,Lfac,Ufac,Pu,nl,nswp,tol,maxiter,stress_data)
penal = 0.5;
% Extract stress data
Bmat = stress_data.B;
Dmat = stress_data.D;
xvec = stress_data.X;
edofMat = stress_data.EDOF;
Emin = stress_data.E(1);
Emax = stress_data.E(2);
% Initialize
r = b - A{1,1}*u;
res0 = norm(b);
% where = stress_data.where;
% where1 = 0;
% maxVM1 = 0;
% maxVM2 = 1000;
% Alternative convergence
dsig1(1:size(xvec),1) = -1000;
dsig2(1:size(xvec),1) = 1000;
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
    % Compute stress and derivative
    STRAIN = u(edofMat)*Bmat';
    STRESS = STRAIN*Dmat;
    VMS2 = STRESS(:,1).^2 - STRESS(:,1).*STRESS(:,2) + STRESS(:,2).^2 + 3*STRESS(:,3).^2;
    VMS = VMS2.^(0.5);
    EVEC = repmat(Emin+(xvec(:).^penal)*(Emax-Emin),1,1);
    TMPGRAD = xvec(:).^(penal-1);
    DEDRHO = penal*(Emax-Emin)*TMPGRAD;
    VM = EVEC.*VMS;
    DSIG = 8*(DEDRHO.*VMS).*(VM.^7);
    altres1 = (DSIG-dsig2)/2;
    altres2 = (DSIG-dsig1);
    scale = max(abs(DSIG));
    altres = max(max(abs(altres1),abs(altres2)))/scale;
%     % Check convergence of max VM and location
%     where2 = where1;
%     where1 = where;
%     [maxVM,where] = max(VM);
%     altres21 = (maxVM-maxVM2)/maxVM/2;
%     altres22 = (maxVM-maxVM1)/maxVM;
%     altres2 = max(altres21,altres22);
%     fprintf(' MGCGiter:%3i relres:%4.3f altres:%4.3f \n',...
%         i,relres,altres);
%     if (where==where1 && where==where2)
%         if max(altres1,altres2) < tol || i>=maxiter
%             break
%         end
%     end
    if altres < tol || i>=maxiter
        break
    end
%     maxVM2 = maxVM1;
%     maxVM1 = maxVM;
    dsig2 = dsig1;
    dsig1 = DSIG;
end
end