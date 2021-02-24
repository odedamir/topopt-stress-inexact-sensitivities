%%%% 2D TOPOLOGY OPTIMIZATION CODE %%%%
% With stress constraints, p-norm approach
% Single density field at eta=0.5 
function minVpnorm_adpt(domain,sizex,sizey,helem,penal,rmin,ft,Sytarget,comptarget,...
    betamax,betamaxsimul,solver,filename)
tic;
close all;
%% MATERIAL PROPERTIES
Emax = 1;
Emin = Emax*1e-6;
nu = 0.3;
%% PREPARE DOMAIN
switch domain
    case 'biclamped'
        [X,T,i_img,j_img] = generate_biclamped(sizex,sizey,helem,false);
    case 'lbracket'
        [X,T,i_img,j_img] = generate_lbracket(sizex,sizey,helem,false);
    case 'ubracket'
        [X,T,i_img,j_img] = generate_ubracket(sizex,sizey,helem,false);
end
%--------------------------------------------------------------
%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nelem = size(T,1);
nodes = size(X,1);
ndof = 2*nodes;
edofMat = zeros(nelem,8);
edofMat(:,1:2) = [2*T(:,1)-1 2*T(:,1)];
edofMat(:,3:4) = [2*T(:,2)-1 2*T(:,2)];
edofMat(:,5:6) = [2*T(:,3)-1 2*T(:,3)];
edofMat(:,7:8) = [2*T(:,4)-1 2*T(:,4)];
iK = reshape(kron(edofMat,ones(8,1))',64*nelem,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelem,1);
% Prologation operators
if solver.type == 4 % Multigrid
    Pu = cell(3,1); % Assumes 4 levels
    Xf = X; % Initialize fine grid
    hgrid = helem;
    for i = 1:3
        [Pu{i,1},Xc] = prepcoarse(Xf,hgrid);
        % Replace fine grid with next level
        Xf = Xc;
        hgrid = hgrid*2;
    end
end
U = zeros(ndof,1); LAM = 0*U;
Bmat = 1/helem*[-1/2 0 1/2 0 1/2 0 -1/2 0
    0 -1/2 0 -1/2 0 1/2 0 1/2
    -1/2 -1/2 -1/2 1/2 1/2 1/2 1/2 -1/2];
% The constitutive matrix
Dmat = 1/(1-nu^2)*[ 1 nu 0;nu 1 0;0 0 (1-nu)/2];
%--------------------------------------------------------------
%% DEFINE LOADS AND SUPPORTS
solids = find(T(:,5)==2); % Find solid elements
padding = find(T(:,6)==1); % Find padding elements
solids_pad = intersect(solids,padding); % Find padding solids
loadedelem = setdiff(solids,solids_pad);
nloadedelem = size(loadedelem,1);
switch domain
    case 'ubracket'
        loadeddof = edofMat(loadedelem,1:2:7);
    otherwise
        loadeddof = edofMat(loadedelem,2:2:8);
end
loadeddof = loadeddof(:);
loadmag = -1; 
load_per_dof = loadmag/nloadedelem/4;
F = sparse(loadeddof,1,load_per_dof,ndof,1);
% Supports
switch domain
    case 'biclamped'
        % Bi-clamped
        [supnodes,~] = find(X(:,1)==0); % Find nodes with x==0
        supdofs1 = union(2*supnodes-1,2*supnodes);
        [supnodes,~] = find(X(:,1)==sizex); % Find nodes with x==sizex
        supdofs2 = 2*supnodes-1;
        supdofs = union(supdofs1,supdofs2);
    case 'lbracket'
        % L-bracket
        [supnodes,~] = find(X(:,2)==sizey); % Find nodes with y==sizey
        supdofs = union(2*supnodes-1,2*supnodes);
    case 'ubracket'
        % U-bracket
        [supnodes,~] = find(X(:,1)==0); % Find nodes with x==0
        supdofs = union(2*supnodes-1,2*supnodes);
end
alldofs = 1:ndof;
freedofs = setdiff(alldofs,supdofs);
% Null space elimination of supports
N = ones(ndof,1); N(supdofs,1) = 0; Null = spdiags(N,0,ndof,ndof);
%--------------------------------------------------------------
%% PREPARE FILTER
iH = ones(nelem*100,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for e1 = 1:nelem
    x1 = T(e1,7);
    y1 = T(e1,8);
    for e2 = 1:nelem
        x2 = T(e2,7);
        y2 = T(e2,8);
        dist = sqrt((x2-x1)^2+(y2-y1)^2);
        if dist<=rmin
            k = k+1;
            iH(k) = e1;
            jH(k) = e2;
            sH(k) = rmin-dist;
        end
    end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
%--------------------------------------------------------------
%% INITIALIZE OPTIMIZATION
maxloop = 90;
volfrac = 0.95;
x = volfrac*ones(nelem,1);
x(T(:,5)==1) = 1e-6; % Voids
x(T(:,5)==2) = 1-1e-6; % Solids
beta = 4;
njumps = 4;
dbeta = (betamaxsimul/beta)^(1/njumps); % Factor for multiplying beta
pace = min(20,maxloop/(njumps+1));
loop = 0;
%% INITIALIZE MMA OPTIMIZER
m     = 2;                % The number of general constraints.
n     = nelem;            % The number of design variables x_j.
xmin  = 1e-6*ones(n,1);   % Column vector with the lower bounds for the variables x_j.
xmin(T(:,5)==2) = 1-1e-3; % Lower bound for solids
xmax  = ones(n,1);        % Column vector with the upper bounds for the variables x_j.
xmax(T(:,5)==1) = 1e-3;   % Upper bound for voids
xold1 = x(:);             % xval, one iteration ago (provided that iter>1).
xold2 = x(:);             % xval, two iterations ago (provided that iter>2).
low   = 0*ones(n,1);      % Column vector with the lower asymptotes from the previous iteration (provided that iter>1).
upp   = ones(n,1);        % Column vector with the upper asymptotes from the previous iteration (provided that iter>1).
a0    = 1;                % The constants a_0 in the term a_0*z.
a     = zeros(m,1);       % Column vector with the constants a_i in the terms a_i*z.
c_MMA = 1000*ones(m,1);   % Column vector with the constants c_i in the terms c_i*y_i.
d     = ones(m,1);        % Column vector with the constants d_i in the terms 0.5*d_i*(y_i)^2.
%% START ITERATION
Sy = Sytarget;
maxVMscaled = 2*Sytarget;
fval = ones(m,1);
Stats = zeros(maxloop,20);
pcgi1 = 0; pcgi2 = 0;
pcgr1 = 0; pcgr2 = 0;
solver.cgtol11 = solver.cgtol1;
solver.cgtol22 = solver.cgtol2;
excessive_pcgi1 = 0;
excessive_pcgi2 = 0;
%% START ITERATION
while ((loop < maxloop || beta<0.99*betamaxsimul) || (maxVMscaled > Sytarget || fval(1,1) > 0))
    loop = loop + 1;
    %% FILTERING OF DESIGN VARIABLES
    if ft == 2
        % Density filter only
        xTilde = (H*x)./Hs;
        xPhys = xTilde;
    elseif ft == 3
        % Density + single projection
        xTilde = (H*x)./Hs;
        xPhys = (tanh(beta*0.5)+tanh(beta*(xTilde-0.5)))/...
            (tanh(beta*0.5)+tanh(beta*(1-0.5)));
    end
    %% FE-ANALYSIS
    sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(Emax-Emin)),64*nelem,1);
    K = sparse(iK,jK,sK);
    K = (K+K')/2;
    % For iterative stress calculation
    stress_data.B = Bmat;
    stress_data.D = Dmat;
    stress_data.X = xPhys;
    stress_data.EDOF = edofMat;
    stress_data.E(1) = Emin;
    stress_data.E(2) = Emax;
    if (solver.type==1)
        % Direct
        U(freedofs) = K(freedofs,freedofs)\F(freedofs);
        flag = 0;
    elseif (solver.type==4)
        % MGCG
        nl = 4;
        Kmg = cell(nl,1);
        Kmg{1,1} = Null'*K*Null - (Null-speye(ndof,ndof));
        for l = 1:nl-1
            Kmg{l+1,1} = Pu{l,1}'*(Kmg{l,1}*Pu{l,1});
        end
        Lfac = chol(Kmg{nl,1},'lower'); Ufac = Lfac';    
        [pcgi1,pcgr1,U] = mgcg_stress(Kmg,F,U,Lfac,Ufac,Pu,nl,1,solver.cgtol11,solver.cgmax1,stress_data);
        flag = 0;
    end
    %%%%%%%%%%%%%%%%
    %% COMPLIANCE AND SENSITIVITY
    ce = sum((U(edofMat)*KE).*U(edofMat),2);
    comp = sum(sum((Emin+xPhys.^penal*(Emax-Emin)).*ce));
    dc = -penal*(Emax-Emin)*xPhys.^(penal-1).*ce;
    dv = ones(nelem,1);
    %% STRESS CONSTRAINT
    % Evaluate stresses
    [appSmax,~,~,~,ADJ,~,dsigvec] = getstress2d(nelem,U,xPhys,edofMat,Emin,Emax,nu,8,helem);
    if (solver.type==1)
        % Direct
        LAM(freedofs) = K(freedofs,freedofs)\ADJ(freedofs);
        pcgi2 = 0; pcgr2 = 0;
    elseif (solver.type==4)
        solver.cgtol22 = solver.cgtol11;
        ADJ = Null*ADJ;
        [pcgi2,pcgr2,LAM] = mgcg_adj(Kmg,ADJ,LAM,Lfac,Ufac,Pu,nl,1,solver.cgtol22,solver.cgmax2,...
            U,edofMat,KE,xPhys,penal,Emax,Emin);
    end
    ce1 = sum(LAM(edofMat)*KE.*U(edofMat),2); % Lam'*K*U
    dsig1 = penal*(Emax-Emin)*xPhys.^(penal-1).*ce1; % Lam'*dKdrho*U
    dsig2 = dsigvec; % dsigaggdrho
    %% FILTERING/MODIFICATION OF SENSITIVITIES
    if  ft == 2
        dv(:) = H*(dv(:)./Hs);
        dsig1(:) = H*(dsig1(:)./Hs);
        dsig2(:) = H*(dsig2(:)./Hs);
        dc(:) = H*(dc(:)./Hs);
    elseif ft == 3
        dxphys = (1 - (tanh(beta*(xTilde(:)-0.5))).^2)*beta / ...
            (tanh(beta*0.5)+tanh(beta*(1-0.5)));
        dv(:) = H*(dv(:).*dxphys(:)./Hs);
        dsig1(:) = H*(dsig1(:).*dxphys(:)./Hs);
        dsig2(:) = H*(dsig2(:).*dxphys(:)./Hs);
        dc(:) = H*(dc(:).*dxphys(:)./Hs);
    end
    dsig = dsig2 + dsig1;
    %% KKT check
    if (loop > 1)
    [~,kktnorm,~] = kktcheck(m,n,xmma,ymma,zmma,lam,...
        xsi,eta,mu_mma,zet,s,xmin,xmax,df0dx,fval,dfdx,a0,a,c_MMA,d);
    else
        kktnorm = 10;
    end
    %% Scaling of constraints (with "dilated" density)
    UxPhys = U;   
    if (1>0)
        % So-called "smooth" scaling
        xDilSharp = (tanh(betamax*0.4)+tanh(betamax*(xTilde-0.4)))/...
            (tanh(betamax*0.4)+tanh(betamax*(1-0.4))); % For measuring "actual" stress
    else
        % So-called "sharp" scaling
        xDilSharp = (tanh(64*0.5)+tanh(64*(xTilde-0.5)))/...
            (tanh(64*0.5)+tanh(64*(1-0.5))); % For measuring "actual" stress
    end
    STRAIN = UxPhys(edofMat)*Bmat';
    STRESS = STRAIN*Dmat;
    Emat = repmat(Emin+xDilSharp(:)*(Emax-Emin),1,3);
    STRESS = STRESS.*Emat;
    VM = STRESS(:,1).^2 - STRESS(:,1).*STRESS(:,2) + STRESS(:,2).^2 + 3*STRESS(:,3).^2;
    VM = VM.^(0.5);
    maxVMscaled = max(VM);
    %% DRAW
    figure(1);
    clf;
    subplot(2,3,1);
    v_img = dsig1;
    top_img = sparse(i_img,j_img,v_img);
    imagesc(top_img);
    axis equal;
    axis tight;
    title('dsig1');
    subplot(2,3,2);
    v_img = dsig2;
    top_img = sparse(i_img,j_img,v_img);
    imagesc(top_img);
    axis equal;
    axis tight;
    title('dsig2');
    subplot(2,3,3);
    v_img = dc;
    top_img = sparse(i_img,j_img,v_img);
    imagesc(top_img);
    axis equal;
    axis tight;
    title('dc');
    subplot(2,3,4);
    v_img = xPhys;
    top_img = sparse(i_img,j_img,v_img);
    imagesc(top_img);
    axis equal;
    axis tight;
    title('xPhys');
    subplot(2,3,5);
    v_img = VM;
    top_img = sparse(i_img,j_img,v_img);
    imagesc(top_img);
    axis equal;
    axis tight;
    title('VM(xDil)');
    drawnow;
    %% METHOD OF MOVING ASYMPTOTES
    xval  = x;
    % minVstCS
    f0val = mean(xPhys);
    if (loop==1)
        scale = 10/f0val;
    end
    df0dx = scale*dv(:)/n;
    if (mod(loop,5) == 1)
        % Scale stress constraint
        ratio = maxVMscaled/Sytarget;
        Sy = appSmax/ratio;
    end
    fval(1,1) = appSmax/Sy - 1;
    dfdx(1,:) = dsig'/Sy;
    if (m > 1) % For testing without compliance
        fval(2,1) = comp/comptarget - 1;
        dfdx(2,:) = dc'/comptarget;
    else
        fval(2,1) = 0;
    end
    fval_mma = fval(1:m,:);
    [xmma,ymma,zmma,lam,xsi,eta,mu_mma,zet,s,low,upp] = ...
        mmasub(m,n,loop,xval,max(xmin,xval-0.1),min(xmax,xval+0.1),xold1,xold2, ...
        f0val,df0dx,fval_mma,dfdx,low,upp,a0,a,c_MMA,d); 
    % Check MMA solution
    if (max(ymma) > 0.1 && solver.cgtol11>2e-6)
        solver.cgtol11 = max(solver.cgtol11/10,1e-6);
        fprintf(' Rejecting approximation, tightening tolerance to %6.3e \n',solver.cgtol11);
        excessive_pcgi1 = excessive_pcgi1 + pcgi1;
        excessive_pcgi2 = excessive_pcgi2 + pcgi2;
        loop = loop-1;
%         % Rollback
%         x = xold1;
%         xold1 = xold2;
%         xold2 = xold3;
%         loop = loop-2;
    else
        if (max(ymma) < 1e-6)
            solver.cgtol11 = min(solver.cgtol11*1.5,1e-1);
            fprintf(' Accepting approximation, loosening tolerance to %6.3e \n',solver.cgtol11);
        else 
            fprintf(' Accepting approximation, keeping tolerance at %6.3e \n',solver.cgtol11);
        end
        % Update MMA Variables
        xnew     = xmma;
        xold3    = xold2(:);
        xold2    = xold1(:);
        xold1    = xval(:);
        change = max(abs(xnew(:)-xval(:)));
        x = xnew;
        %% CONTINUATION
        if (mod(loop,pace) == 0)
            beta = min(beta*dbeta,betamaxsimul);
        end
        
        %% PRINT RESULTS
        fprintf(' It.:%3i Vol.:%4.3f ch.:%4.3f VMapp:%6.3e VMacc:%6.3e g: %6.3e %6.3e betaHS:%6.3f\n',...
            loop,mean(xPhys(:)),change,appSmax,maxVMscaled,fval',beta);
        fprintf(' State solver %4i %6.3e Adjoint solver %4i %6.3e\n',...
            pcgi1,pcgr1,pcgi2,pcgr2);
        %% KEEP STATS
        Stats(loop,1:15+2*(m-1)+1) = [f0val fval' appSmax maxVMscaled beta pcgi1 pcgr1 pcgi2 pcgr2...
            ymma' zmma' lam' flag kktnorm solver.cgtol11]; 
    end   
end
runtime = toc;
save(filename);
end



