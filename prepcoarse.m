%% FUNCTION prepcoarse - PREPARE MG PROLONGATION OPERATOR
function [Pu,Xc] = prepcoarse(Xf,h)
% Create coarse grid
Xc = 0*Xf;
count = 0;
for i = 1:size(Xf,1)
    if(mod(Xf(i,1),2*h) == 0)
        if(mod(Xf(i,2),2*h) == 0)
            % This node is on the coarse mesh
            count = count + 1;
            Xc(count,1:2) = Xf(i,1:2);
        end
    end
end
Xc(count+1:end,:) = [];
% Create prolongation operator
nnodesf = size(Xf,1);
nnodesc = size(Xc,1);
maxnum = nnodesf*20;
iP = zeros(maxnum,1); jP = zeros(maxnum,1); sP = zeros(maxnum,1);
% Weights for fixed distances to neighbors on a structured grid 
vals = [1,0.5,0.25];
count = 0;
for nc = 1:nnodesc
    col = nc; % The columns of P are related to coarse node
    for nf = 1:nnodesf
        row = nf; % The rows of P are related to fine node
        % Find distance between nodes
        dist = sqrt((Xc(nc,1)-Xf(nf,1))^2 + (Xc(nc,2)-Xf(nf,2))^2);
        if (dist <= h*1.5) % Is within influence
            ind = ceil(1 + dist/h);
            count = count+1; iP(count) = 2*row-1; jP(count) = 2*col-1; sP(count) = vals(ind);
            count = count+1; iP(count) = 2*row; jP(count) = 2*col; sP(count) = vals(ind);
        end
    end
end
% Assemble matrices
Pu = sparse(iP(1:count),jP(1:count),sP(1:count));
end