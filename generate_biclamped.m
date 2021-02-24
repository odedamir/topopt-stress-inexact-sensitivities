function [X,T,row,col] = generate_biclamped(sizex,sizey,helem,doplot)
% Define grid parameters
dx = helem;
dy = helem;
pad = 0;
% Create nodal grid for FEA and optimization
[Xnode,Ynode] = meshgrid(-pad:dx:sizex+pad,-pad:dy:sizey);
nodelist(:,1) = Xnode(:);
nodelist(:,2) = Ynode(:);
xnodeout = [];
ynodeout = [];
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
[Xelem,Yelem] = meshgrid(-pad+dx/2:dx:sizex+pad-dx/2,-pad+dy/2:dy:sizey-dy/2);
elemlist(:,1) = Xelem(:);
elemlist(:,2) = Yelem(:);
xelemout = [];
yelemout = []; 
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
    % Assign solids (for load)
    if (x_cent>sizex-5 && y_cent< 5); T(e,5:6) = [2 0]; end
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
xmin = -pad+dx/2; 
ymin = -pad+dy/2; ymax = sizey-dy/2;
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