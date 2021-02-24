%% FUNCTION VCycle - COARSE GRID CORRECTION
function z = VCycle(A,r,Lfac,Ufac,Pu,l,nl,invD,omega,nswp)
z = 0*r;
z = smthdmpjac(z,A{l,1},r,invD{l,1},omega,nswp);
Az = A{l,1}*z;
d = r - Az;
dh2 = Pu{l,1}'*d;
if (nl == l+1)
    vh2 = Ufac \ (Lfac \ dh2);
else
    vh2 = VCycle(A,dh2,Lfac,Ufac,Pu,l+1,nl,invD,omega,nswp);
end
v = Pu{l,1}*vh2;
z = z + v;
z = smthdmpjac(z,A{l,1},r,invD{l,1},omega,nswp);
end