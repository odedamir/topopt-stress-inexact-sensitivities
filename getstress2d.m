function [appSmax,trueSmax,where,SumappSmax,ADJ,vmvec,dsigvec] = ...
    getstress2d(n,U,x,edofMat,Emin,Emax,nu,PN,side)
penal = 0.5;
% The B-matrix
Bmat = 1/side*[-1/2 0 1/2 0 1/2 0 -1/2 0
    0 -1/2 0 -1/2 0 1/2 0 1/2
    -1/2 -1/2 -1/2 1/2 1/2 1/2 1/2 -1/2];
% The constitutive matrix
D = 1/(1-nu^2)*[ 1 nu 0;nu 1 0;0 0 (1-nu)/2];
% The von Mises Vmatrix
V = [1 -0.5 0;
    -0.5 1  0;
    0 0 3];
% Loop over all elements
SumappSmax = 0;
ADJ = 0*U;
vmvec = 0*x;
dsigvec = vmvec;
for i = 1:n
    edof = edofMat(i,:);
    Ue = U(edof,1);
    strain = Bmat*Ue;
    stress = D*strain; % Assuming E = 1
    vms2 = stress'*V*stress; % This is von Mises stress ^2
    vms = sqrt(vms2); % This is the von Mises stress for E=1
    Ee = Emin+(Emax-Emin)*x(i,1)^penal; % Actual E at 
    tmpgrad = x(i,1)^(penal-1);
    dEedrho = penal*(Emax-Emin)*tmpgrad;
    vm = Ee*vms; % This is the macroscopic SIMP von Mises stress
    vmvec(i,1) = vm; % Save in vector
    % Absolute measure
    SumappSmax = SumappSmax + vm^PN; % Collect relative contributions in P-norm
    % Compute contribution to adjoint
    SVCB = stress'*V*D*Bmat; % [CBU]'*[V]*[CB]
    adje = Ee*1/vms*SVCB;
    adje = adje*PN*vm^(PN-1);
    dsige = dEedrho*vms*PN*vm^(PN-1);
    dsigvec(i,1) = dsige;
    ADJ(edof,1) = ADJ(edof,1) - adje';
end
% Scale ADJ
ADJ = ADJ*1/PN*(SumappSmax^(1/PN-1));
dsigvec = dsigvec*1/PN*(SumappSmax^(1/PN-1));
appSmax = SumappSmax^(1/PN);
[trueSmax,where] = max(vmvec);
