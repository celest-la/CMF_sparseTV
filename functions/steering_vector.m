function steer_vec=steering_vector(X,Z,h,f0)

    dist=sqrt((X(:)-h.xp).^2 + Z(:).^2);
    tau_pix = dist./h.c0;  
    steer_vec=1/h.Ne*exp(-2*1i*pi*f0*tau_pix');

