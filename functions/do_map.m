function [map] = do_map(algo, CSM, axf, X, Z, h, varargin)

%-- INPUTS
%   algo : name of the algorithme to compute
%       'pci'  : FD-DAS
%       'rcb'  : Robust Capon Beamforming, takin parameter 'eps' (optional input)
%   CSM : Cross-Spectral Matrix [elem elem]
%   axf : frequency bins
%   X,Z : Grid of reconstruction area, from meshgrid function
%   h : probe header
%       h.xp : position of the element on the x axis (m)
%       h.c0 : speed of sound (m/s)
%
%-- (optional) INPUTS
%   eps : (default to 1) parameter for RCB method
%   noWb : (default to 1) to display or not waitbar during RCB computation
%       0 : display
%       1 : no display 

%%

p = inputParser;

addRequired(p,'algo');
addRequired(p,'CSM');
addRequired(p,'axf');
addRequired(p,'X');
addRequired(p,'Z');
addRequired(p,'h');

addParameter(p,'eps',1);
addParameter(p,'noWb',1);

parse(p,algo,CSM,axf,X,Z,h,varargin{:});

name_algo = p.Results.algo;
noWb = p.Results.noWb;

%% Compute distance and delays for Steering Vector
% Number of sensors 
M = size(CSM,1);

% Compute the distance from each pixel to each elements
d_pix = sqrt((X(:)-h.xp).^2 + Z(:).^2);     % [Pixels, Elements]

% Compute the corresponding delay
tau_pix = d_pix./h.c0;

%% Reconstruction Part

% Init. cavitation maps
map = zeros([size(X),length(axf)]);
map_t = zeros([length(tau_pix),1]); %Temporary Map
%Lambda = zeros([size(X),length(axf)]);

% Init. weight vectors and other values for each algorithm
switch name_algo
    case 'rcb'
        w   = zeros([size(X), M]);
        I = eye(M);  % Identity Matrix
    case 'pci'
        w = 1/M;  %[Elements,1]

end

%-- loop on frequencies 
for nf = 1:length(axf)
    
    nCSM = squeeze(CSM(:,:,nf));
    
    switch name_algo
        
        case 'pci'
%===========================================================
            W = nCSM;
            
            %Steering vector (at pixel of interest)
            ab = w.*exp(-2*1i*pi*axf(nf)*tau_pix');
            
            %Normalization
            n_ab = vecnorm(ab,2,1);
            ab   = ab./ n_ab;
            
            %Compute Power
            map_t = real( sum((ab'*W) .* (ab.'),2) ./ (n_ab').^2);
            
            %Reshape the Power Matrix map_t
            map(:,:,nf)     = reshape(map_t ,size(X));
            
             
        case 'rcb'
%===========================================================
            eps = p.Results.eps;
            tau_pix = reshape(tau_pix,[size(X),M]);
            % Compute the eigenvalue decomposition of CSM
            
            if cond(nCSM) == 0
                error("La matrice est non inversible car det(CSM)=0");
            end
            
            [U,V]= eig(nCSM);  
            if noWb == 0; wb = waitbar(0,"rcb en cours"); end
            for ix = 1:size(X,2)
                for iz = 1:size(X,1)
                    
                    % Select the steering vector for the Pixel of Interest
                    ab = exp(-2*1i*pi*axf(nf)*squeeze(tau_pix(iz,ix,:)));
                    
                    if eps~=0
                        %- Vectors used later 
                        sumU = U'*(ab);
                        valp = diag(V);
                        
                        %- Theoretical Boundaries for Lambda
                        % The Low One
                        lb = (norm(ab) - sqrt(eps)) / (max(valp)*sqrt(eps));

                        % Find the optimal Lambda: Newton-Raphson method
                        % Initialization lambda and number of iteration
                        lambda = 0.5*lb;
                        itr = 0;

                        % Find the optimal lambda in few iterations
                        while (itr < 20)
                            % Fonction f(lambda)=0 to be solved
                            f = sum(abs(sumU).^2 ./ (1 + lambda * valp).^2)-eps;
                            % Derivative function of f(lambda)
                            fp = -2 * sum((valp .* abs(sumU).^2) ./ (1 + lambda * valp).^3);
                            % Update lambda, and number of iteration
                            lambda = lambda - f / fp;
                            itr = itr + 1;
                            
                        end
                   
                     map(iz,ix,nf) = 1 / ((ab')*U*V* ...
                        ((lambda^(-2)*I+(2/lambda)*V+V*V)\I)*(U')*ab);
                    

    
                    %- Compute the norm of ac
                    %  the steering vector comprising uncertainities 
                    ac = ab - (U*((I+lambda*V)\I)*U')*ab;  

                    %- Correct scaling ambiguity on Power estimates with RCB
                    map(iz,ix,nf) = map(iz,ix,nf)* norm(ac).^2 / M; 

                    else
                        iCSM    = inv(nCSM);
                        map(iz,ix,nf) = 1 / (ab' * iCSM * ab);
                    end
                    if noWb == 0; waitbar(((ix-1)*size(X,1) + iz)/(size(X,1)*size(X,2)),wb,"rcb en cours"); end
                end
            end
            if noWb == 0; close(wb); end
    end
end
      
map = real(map);
%End of function
end
        