function [X,T,row,col] = generate_ubracket(sizex,sizey,helem,doplot)
% Define grid parameters
dx = helem;
dy = helem;
pad = 8; 
% Create nodal grid for FEA and optimization
[Xnode,Ynode] = meshgrid(0:dx:sizex,-pad:dy:sizey+pad);
nodelist(:,1) = Xnode(:);
nodelist(:,2) = Ynode(:);
xnodeout1 = find(nodelist(:,1)>0.5*sizex+pad);
xnodeout2 = find(nodelist(:,1)<0.75*sizex-pad);
xnodeout = intersect(xnodeout1,xnodeout2);
ynodeout = find(nodelist(:,2)>0.5*sizey+pad);
xynodeout = intersect(xnodeout,ynodeout);
nodelist_clean = nodelist;
nodelist_clean(xynodeout,:) = [];

% Plot nodes as points
if (doplot)
    plot(nodelist_clean(:,1),nodelist_clean(:,2),'o');
    hold on;
    axis equal
    axis tight
end

% Create element grid for FEA and optimization
[Xelem,Yelem] = meshgrid(dx/2:dx:sizex-dx/2,-pad+dy/2:dy:sizey+pad-dy/2);
elemlist(:,1) = Xelem(:);
elemlist(:,2) = Yelem(:);
xelemout1 = find(elemlist(:,1)>0.5*sizex+pad);
xelemout2 = find(elemlist(:,1)<0.75*sizex-pad);
xelemout = intersect(xelemout1,xelemout2);
yelemout = find(elemlist(:,2)>0.5*sizey+pad);
xyelemout = intersect(xelemout,yelemout);
elemlist_clean = elemlist;
elemlist_clean(xyelemout,:) = [];

% Create element connectivity
nelem = size(elemlist_clean,1);
T = zeros(nelem,8);
for e = 1:nelem
    % Centroid coordinates
    x_cent = elemlist_clean(e,1);
    y_cent = elemlist_clean(e,2);
    T(e,7:8) = [x_cent y_cent];
    % Nodes
    x1 = find(nodelist_clean(:,1)==(x_cent-dx/2));
    y1 = find(nodelist_clean(:,2)==(y_cent-dy/2));
    n1 = intersect(x1,y1);
    x2 = find(nodelist_clean(:,1)==(x_cent+dx/2));
    y2 = find(nodelist_clean(:,2)==(y_cent-dy/2));
    n2 = intersect(x2,y2);
    x3 = find(nodelist_clean(:,1)==(x_cent+dx/2));
    y3 = find(nodelist_clean(:,2)==(y_cent+dy/2));
    n3 = intersect(x3,y3);
    x4 = find(nodelist_clean(:,1)==(x_cent-dx/2));
    y4 = find(nodelist_clean(:,2)==(y_cent+dy/2));
    n4 = intersect(x4,y4);
    T(e,1:4) = [n1 n2 n3 n4];
    % Assign voids
    if (y_cent<0); T(e,5:6) = [1 1]; end
    if (y_cent>sizey); T(e,5:6) = [1 1]; end
    if (x_cent>0.5*sizex && x_cent<0.75*sizex && y_cent>0.5*sizey); T(e,5:6) = [1 1]; end
    % Assign solids (for load)
    if (x_cent>sizex-5 && y_cent>sizey-5); T(e,5:6) = [2 0]; end
    if (x_cent>sizex-5 && y_cent>sizey); T(e,5:6) = [2 1]; end
end

% Plot elements
if (doplot)
    for e = 1:nelem
        if (T(e,5) == 0) % Regular element
            plot(T(e,7),T(e,8),'ro'); 
        end
        if (T(e,5) == 2) % Solid element
            plot(T(e,7),T(e,8),'k*'); 
        end
        if (T(e,5) == 1) % Void element
            plot(T(e,7),T(e,8),'m+'); 
        end
    end
end

% Create matrix representation of topology
xmin = dx/2; 
ymin = -pad+dy/2; ymax = sizey+pad-dy/2;
nrows = (ymax-ymin)/dy + 1;
col = zeros(nelem,1);
row = col;
val = col;
for e = 1:nelem
    col(e) = (T(e,7)-xmin)/dx + 1;
    row(e) = nrows - ((T(e,8)-ymin)/dy);
    val(e) = T(e,5);
end

% For output
X = nodelist_clean;