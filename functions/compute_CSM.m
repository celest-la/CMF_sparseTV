function [CSM,axf,S] = compute_CSM(raw,fs,K,overlap,f,varargin)
% Function computing the Discrete Fourier Transform for K overlapping
% snapshot of the raw signals and the CSM matrix
%-- INPUTS:
%   raw : matrix of N raw signals recorded by the array [Time, Elements]
%   fs  : frequency sampling of rawdata
%   K   : number of snapshots, typically K >= number of elements
%   overlap: overlap percentage of snapshot length (%) typically 0.9
%   f : frequencies on which recontruction is made
%
%-- (optional) INPUTS :
%   Nzp : (default, 0) Number of samples on which the fft is computed
%           if 0, the value is taken to have a frequency resolution of 1%
%           of f.
%   nBin : (default 1) Number of frequency to take around the frequency f
%           (not tested for length(f)>1).
%   doDisp : (default "on") To display or not the warnings if Tsnaps < 20
%           periods of f, the duration of Tsnaps and the number of periods of f in
%           Tsnaps.
%
%-- OUTPUT :
%   CSM is a [Ne x Ne x Nf] matrix, containing the CSM at each frequency
%   bin. Ne is the number of sensors, Nf is the number of frequency bins
%   in w.
%   axf is the value of the frequency bin.
%
       


%% 
p = inputParser;

addRequired(p,'raw');
addRequired(p,'fs');
addRequired(p,'K');
addRequired(p,'overlap');
addRequired(p,'f');

addParameter(p,'Nzp',0);
addParameter(p,'nBin',1);
addParameter(p,'doDisp',"on");

parse(p,raw,fs,K,overlap,f,varargin{:});

nBin = p.Results.nBin;

%% Function
%-----------------------------------------------------------

% Define the variables
    [Nt, Ne] = size(raw);                                      
    % Nt : Time samples
    % Ne : Number of elements
    nw = floor(Nt/( K*(1-overlap)+overlap));
    % nw : number of samples in a snapshot
   
    lag = round((1-overlap)*nw);         
    % lag : number of samples between snapshots
    if nw == Nt % Only 1 snapshot
        start = 1;
    else
        start = 1:lag:Nt-nw;                   
    end
    % start : vector containing the starting samples of each snapshot
    nSnapshots = length(start);            
    % nSnapshots : number of snapshots, normally equivalent to K
    
    
    %==========
    % Check the duration of a snapshot according to the reconstruction
    % frequency f, suppposed to contain at least 20 periods of f.
    % Tsnaps : duration of a snapshot
    Tsnaps = nw/fs; 
    if Tsnaps*f <= 20
        if p.Results.doDisp == "on"
            warning(['Tsnaps (' num2str(Tsnaps) ') est peut-être trop court pour avoir une bonne reconstruction. Il y a moins de 20 périodes dans un snapshots. Choisir K plus petit, ou overlap plus grand.']);
        end
    end
    if p.Results.doDisp == "on"
        disp([' Tsnaps = ' num2str(Tsnaps*1e3) ' ms soit ' num2str(Tsnaps*f) ' périodes.'])
    end
    
    %==========
    % Check if Shannon is respected
    if fs < 2*f % Shannon's condition
        error("fs < 2f : La condition de Shannon n'est pas respectée.")
    end

    %==========
    %  Set Nzp
    Nzp = p.Results.Nzp;
    if Nzp == 0 % If Nzp is not set by the user
        deltaf = 0.01*min(f);       % Frequency resolution of fft set as 1% of f
        Nzp = round(fs/deltaf);     % Number of samples on which fft is computed to insure deltaf
        if Nzp<nw; Nzp = nw; end    % Insure that Nzp is superior to nw
    end
    Nzp = 2^nextpow2(Nzp);          % Set Nzp as a power of 2 value.
 

    %==========
    % Find indices of frequencies in frequency vector
    vec_freq = ((1:Nzp)-1)/Nzp * fs;                    % Frequency vector
    ind = zeros([length(f) nBin]);                      % Initialize indice matrix [number of f, number of Bin for each f]
    for i = 1:length(f) 
        [~,ii]= min(abs(vec_freq-f(i)));                % Find indice correspond to f
        if nBin == 1; ind(i,:) = ii;                    % If 1 bin per frequency, fill ind
        else; ind(i,:) = ii-nBin/2:ii+nBin/2-1; end     % If more than 1 bin, take nBin around f
        axf(i) = vec_freq(ind(i,:));                    % Vector containing the frequency bins.
    end
    
    %==========
    %-- Create a window to smooth edges of snapshots (hann window)
    win = repmat(tukeywin(nw,1),[1,Ne]); 
     
    %==========
    % Truncate rawdata and compute the fft.

    TF_f = zeros([Ne,length(ind),nSnapshots]); 
    % TF_f : matrix containing the fft of the frequency bins of each snapshots
    for i = 1:nSnapshots
        %-- Select K blocs of overlapping truncated signals
        raw_snap = raw(start(i) + (0:nw-1),:).*win;
        %-- Compute the DFT of each Snapshot
 
        %Original method
      TF = fft(raw_snap,Nzp); 
        %-- Take the values of the frequency bins only
        TF_f(:,:,i) = TF(ind,:).'; 

    end
    %==========  
    % Compute the CSM matrix
    CSM     = zeros(Ne,Ne,length(ind)); % Averaged CSM
    % CSM : total CSM matrix, average of all CSM from all snapshots
    for i=1:length(ind)
        S = squeeze(TF_f(:,i,:));
        CSM(:,:,i) = S*S'/nSnapshots;
    end
        


end

