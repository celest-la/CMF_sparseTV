%%%ADMM Function built from Boyd et al. 2010
%%%adapted for generalized lasso with the addition of the L2 regularization
%%%term
function [x,history,inv_op]=generalizedLassoAdmm(y,A,lambda,lambdal2,rho,F,inv_op)
%%%% Initialisation
    F=sparse(F);
    max_iter=1000;
    n=size(F,1);
    if numel(inv_op)==0
        inv_op=inv((A'*A)+rho*(F'*F)+lambdal2*eye(size(A,2)));
    end

    if canUseGPU() 
        F=gpuArray(F);
        A=gpuArray(A);
        y=gpuArray(y);
        lambdal2=gpuArray(lambdal2);
        inv_op=gpuArray(inv_op);
        z=zeros(n,1,'gpuArray');
        u=zeros(n,1,'gpuArray');
        x=zeros(size(F,2),1,'gpuArray');    
        alpha=gpuArray(1.5);
        lambda=gpuArray(lambda);
        ABSTOL   = gpuArray(1e-6);
        RELTOL   = gpuArray(1e-2);
    
    else
        z=zeros(n,1);
        u=zeros(n,1);
        x=zeros(size(F,2),1);
        alpha=1.5;
         ABSTOL   = 1e-6;
        RELTOL   = 1e-2;
    end
    
    Aty=A'*y;
    
    for k=1:max_iter
        x=inv_op*(Aty+rho*F'*(z-u));
        zold=z;
        Ax_hat=alpha*F*x+(1-alpha)*zold;
        z=shrinkage(Ax_hat+u,lambda/rho);
        u=u+(Ax_hat-z);
        
        %%%Stopping criteria
        history.objval(k)  = objective(y, lambda,lambdal2, A, x, z);
        
        history.r_norm(k)  = norm(F*x - z);
        history.s_norm(k)  = norm(-rho*F'*(z - zold));
        history.eps_pri(k) = sqrt(n)*ABSTOL + RELTOL*max(norm(F*x), norm(-z));
        history.eps_dual(k)= sqrt(n)*ABSTOL + RELTOL*norm(rho*F'*u);

            if (history.r_norm(k) < history.eps_pri(k) && history.s_norm(k) < history.eps_dual(k))
                 break;
            end
            gather(x);
            gather(history);
        
    end 
end

function obj = objective(b, lambda,lambdal2,A, x, z)
    obj = 0.5*norm(A*x - b,2) +0.5*lambdal2*norm(x,2) +norm(lambda.*z,1);
end

function y = shrinkage(a, kappa)
    y = max(0, a-kappa) - max(0, -a-kappa);
end
