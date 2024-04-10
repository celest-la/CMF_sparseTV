%%%%%%%%%%%%% CMF%%%%
%%This function is inspired from the work of Sarradj (2012) on Cross-Matrix
%%fitting. 
%%RF : the raw data, 
% f0 : the imaging frequency, h : the information about
%%the probe and the medium, X and Z the grid, K and overlap the snapshot

function [map,history,inv_op,A,C,np]=CMF(A,C,X,lambdaTV_x,lambdaTV_z,lambdaL1,lambdaL2,rho,inv_op)
    
    
    A=normalize(A,'norm');
    C=normalize(C,'norm');
    
    max_lambda=norm(A'*C,100000); 
    
    m=size(A,2);
    
    %%% Operator for the Vertical total variation
    dv=[1 -1];
    DH=convmtx2(dv,size(X,1),size(X,2));
    
    %%% Operator for the horizontal total variation
    dh=[ 1 ;-1];
    DV=convmtx2(dh,size(X,1),size(X,2));
    
    
    %%% Operator for the sparsity
    D1=eye(m);
    
    
    D=[DH;DV;D1];
    
    lambda_gen=[repmat(lambdaTV_x*max_lambda,size(DH,1),1);repmat(lambdaTV_z*max_lambda,size(DV,1),1);repmat(max_lambda*lambdaL1,size(D1,1),1)];
    tic
    [alpha_gpu,history_gpu,inv_op_gpu]=generalizedLassoAdmm(C,A,lambda_gen,lambdaL2*max_lambda,rho,D,inv_op);
    toc
    alpha=gather(alpha_gpu);
    history=gather(history_gpu);
    inv_op=gather(inv_op_gpu);
    map=reshape(alpha,size(X));
    map=full(map);

