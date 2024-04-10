function [A,C,axf]=getAC(RF,fim,h,X,Z,K,overlap)
%%% The function return the CSM reorganized and create the steering
%%% operator 
%%RF: RF signals to beamform
%%%fim: frequence of imaging
%%%h: parameters of the probe
%%% X,Z, grid of the pixels
%%% K: number of snapshots
%%% overlap: overlapping of the snaphots
np=length(X(:)); 

%% CSM and steering vector acquisition

[CSM,axf] = compute_CSM(RF,h.fs,K,overlap,fim);
steer_vec=steering_vector(X,Z,h,axf);

%% Shaping of the CSM and A for inverse problem :

steer_conj=steer_vec';
steer_conj=reshape(steer_conj,1,np,h.Ne);
steer_vec=reshape(steer_vec,h.Ne,np,1);
Bc=steer_vec.*steer_conj;
Bctemp=permute(Bc,[1 3 2]);
Ac=reshape(Bctemp,h.Ne*h.Ne,np); %%H.Ne: number of elements, np : size of the grid

%%We keep only the upper triangle, because the matrix is symmetric. and diagonal removal 

%Triangle Ac:
mask=triu(ones(h.Ne,h.Ne),0)>0;
mask=reshape(mask',h.Ne*h.Ne,1);
A=complexToReal(Ac(mask,:));


%Triangle CSM :
CSMr=CSM(:);
C=complexToReal(CSMr(mask,:));
end

function rX=complexToReal(X)

X1=real(X);
X2=imag(X);
rX=[X1;X2];
end

function steer_vec=steering_vector(X,Z,h,f0)

    dist=sqrt((X(:)-h.xp).^2 + Z(:).^2);
    tau_pix = dist./h.c0;  
    steer_vec=1/h.Ne*exp(-2*1i*pi*f0*tau_pix');

end

%%This function separate the real and im part of a matrice

