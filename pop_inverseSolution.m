function EEG = pop_inverseSolution(EEG, windowSize, overlaping, solverType, saveFull, account4artifacts, postprocCallback)
if nargin < 1, error('Not enough input arguments.');end
if nargin < 5
    answer = inputdlg({'Window size','Overlaping (%)', 'Save full PCD', 'Solver type','Account for artifacts'},'pop_inverseSolution',1,{num2str((40/1000)*EEG.srate),'50', 'bsbl', 'true', 'true'});
    if isempty(answer)
        error('Not enough input arguments.');
    else
        windowSize = str2double(answer{1});
        overlaping = str2double(answer{2});
        solverType = lower(answer{3});
        saveFull = str2num(lower(answer{4})); %#ok
        account4artifacts = str2num(lower(answer{5})); %#ok
    end
end
if ~isnumeric(windowSize)
    disp('Invalid input for windowSize parameter, we will use the default value.')
    windowSize= (40/1000)*EEG.srate;
end
if ~isnumeric(overlaping)
    disp('Invalid input for overlaping parameter, we will use the default value.')
    overlaping= 50;
end
if ~any(ismember({'loreta','bsbl'},solverType))
    disp(['Solver ' solverType 'is not available, we will use bsbl.']);
    solverType = 'bsbl';
end
if ~islogical(saveFull)
    disp('Invalid input for saveFull parameter, we will use the default value.')
    saveFull= true;
end
if ~islogical(saveFull)
    disp('Invalid input for account4artifacts parameter, we will use the default value.')
    account4artifacts= true;
end
if nargin < 7, postprocCallback = [];end
overlaping = overlaping/100;

% Load the head model
try
    hm = headModel.loadFromFile(EEG.etc.src.hmfile);
catch
    h = errordlg('EEG.etc.src.hmfile seems to be corrupted or missing, to set it right next we will run >> EEG = pop_forwardModel(EEG)');
    waitfor(h);
    EEG = pop_forwardModel(EEG);
    try
        hm = headModel.loadFromFile(EEG.etc.src.hmfile);
    catch
        errordlg('For the second time EEG.etc.src.hmfile seems to be corrupted or missing, try the command >> EEG = pop_forwardModel(EEG);');
        return;
    end
end

% Select channels
labels_eeg = {EEG.chanlocs.labels};
[~,loc] = intersect(lower(labels_eeg), lower(hm.labels),'stable');
EEG = pop_select(EEG,'channel',loc);

% Initialize the inverse solver
Ndipoles = size(hm.cortex.vertices,1);
if account4artifacts && exist('Artifact_dictionary.mat','file')
    [H, Delta, blocks, indG, indV] = buildAugmentedLeadField(hm);
else
    norm_K = norm(hm.K);
    H = hm.K/norm_K;
    Delta = hm.L/norm_K;
    H = bsxfun(@rdivide,H,sqrt(sum(H.^2)));
    blocks = hm.indices4Structure(hm.atlas.label);
    indG = (1:size(H,2))';
    indV = [];
end
Nx = size(H,2);
solver = PEB(H, Delta, blocks);
options = solver.defaultOptions;
options.verbose = false;
if strcmpi(solverType,'loreta')
    options.doPruning = false;
end
EEG.data = double(EEG.data);
Nroi = length(hm.atlas.label);
try
    X = zeros(Nx, EEG.pnts, EEG.trials);
catch ME
    disp(ME.message)
    disp('Using a LargeTensor object...')
    try
        X = LargeTensor([Nx, EEG.pnts, EEG.trials]);
    catch
        disp('Not enough disk space to save src data on your tmp/ directory, we will try your home/ instead.');
        [~,fname] = fileparts(tempname);
        if ispc
            homeDir = getenv('USERPROFILE');
        else
            homeDir = getenv('HOME');
        end
        filename = fullfile(homeDir,fname);
        X = invSol.LargeTensor([Nx, EEG.pnts, EEG.trials], filename);
    end
end
X_roi = zeros(Nroi, EEG.pnts, EEG.trials);

% Construct the average ROI operator
P = hm.indices4Structure(hm.atlas.label);
P = double(P);
P = sparse(bsxfun(@rdivide,P, sum(P)))';

% Check if we need to integrate over Jx, Jy, Jz components
if Nx == Ndipoles*3
    P = [P P P];
    isVect = true;
else
    isVect = false;
end

prc_5 = round(linspace(1,EEG.pnts,30));
iterations = 1:windowSize:EEG.pnts-windowSize;
prc_10 = iterations(round(linspace(1,length(iterations),10)));

windowSize=max([5,windowSize]);
smoothing = filterDesign(EEG.srate);

logE = zeros([length(1:windowSize:EEG.pnts),EEG.trials]);
gamma = zeros([solver.Ng,length(1:windowSize:EEG.pnts),EEG.trials]);
indGamma = EEG.times(1:windowSize:EEG.pnts);

% Perform source estimation
fprintf('PEB source estimation...\n');

for trial=1:EEG.trials
    fprintf('Processing trial %i of %i...',trial, EEG.trials);
    c = 1;
    for k=1:windowSize:EEG.pnts
        loc = k:k+windowSize-1;
        loc(loc>EEG.pnts) = [];
        if isempty(loc), break;end
        if length(loc) < windowSize
            [X(:,loc(1):end,trial),~,~,gamma(:,c,trial), logE(c,trial)] = solver.update(EEG.data(:,loc(1):end,trial), [],[],options);
            X_roi(:,loc(1):end,trial) = computeSourceROI(X(indG,loc(1):end,trial),P,isVect);
            break;
        end
        
        % Source estimation
        [X(:,loc,trial),~,~,gamma(:,c,trial), logE(c,trial)] = solver.update(EEG.data(:,loc,trial),[],[],options);
        indGamma(c) = loc(end);
        
        % Compute average ROI time series
        X_roi(:,loc,trial) = computeSourceROI(X(indG,loc,trial),P,isVect);
        
        % Post-processing (if any)
        if ~isempty(postprocCallback)
            EEG = postprocCallback(EEG, Gamma, EEG.times(loc(end)), trial);
        end
        % Progress indicatior
        [~,ind] = intersect(loc(1:windowSize),prc_5);
        if ~isempty(ind), fprintf('.');end
        prc = find(prc_10==k);
        if ~isempty(prc), fprintf('%i%%',prc*10);end
        c = c+1;
    end
    X(:,:,trial) = filtfilt(smoothing,1,X(:,:,trial)')';
    X_roi(:,:,trial) = filtfilt(smoothing,1,X_roi(:,:,trial)')';
    fprintf('\n');
end
fprintf('done\n');
EEG.etc.src.act = X_roi;
EEG.etc.src.roi = hm.atlas.label;
EEG.etc.src.gamma = gamma;
EEG.etc.src.indGamma = indGamma;
EEG.etc.src.H = H;
EEG.etc.src.indG = indG;
EEG.etc.src.indV = indV;
EEG.etc.src.logE = logE;

fprintf('Cleaning data...');
dim = size(X);
if length(dim) < 3
    EEG.data = H(:,indG)*X(indG,:);
else
    EEG.data = reshape(H(:,indG)*reshape(X(indG,:,:),length(indG),[]),[EEG.nbchan, EEG.pnts, EEG.trials]);
end
fprintf('done.\n');

if saveFull
    try
        EEG.etc.src.actFull = X;
    catch
        EEG.etc.src.actFull = invSol.LargeTensor([Nx, EEG.pnts, EEG.trials], tempname);
    end
else
    EEG.etc.src.actFull = [];
end
EEG.history = char(EEG.history,['EEG = pop_inverseSolution(EEG, ' num2str(windowSize) ', ' num2str(overlaping) ,', ' solverType ', ' num2str(saveFull) ', ' num2str(account4artifacts) ');']);
disp('The source estimates were saved in EEG.etc.src');
end

%%
function x_roi = computeSourceROI(x,P,isVect)
if isVect
    x_roi = sqrt(P*(x.^2));
else
    x_roi = P*x;
end
end

%%
function b = filterDesign(Fs,N)
if nargin < 2
    N = 16;              % Order;
end
    
% FIR Window Lowpass filter designed using the FIR1 function.
% All frequency values are in Hz.

Fc   = floor(Fs/2)-20;  % Cutoff Frequency
flag = 'scale';         % Sampling Flag
Beta = 0.5;             % Window Parameter

% Create the window vector for the design algorithm.
win = kaiser(N+1, Beta);

% Calculate the coefficients using the FIR1 function.
b  = fir1(N, Fc/(Fs/2), 'low', win, flag);
end
% [EOF]
