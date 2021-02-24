%% FUNCTIODN smthdmpjac - DAMPED JACOBI SMOOTHER
function [u] = smthdmpjac(u,A,b,invD,omega,nswp)
for i = 1:nswp
    u = u - omega*invD.*(A*u) + omega*invD.*b;
end
end