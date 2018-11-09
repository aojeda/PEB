function EEG = pop_inverseSolution(EEG, windowSize, overlaping, solverType, saveFull, account4artifacts, postprocCallback)
persistent solver

if nargin < 1, error('Not enough input arguments.');end
if nargin < 5
    answer = inputdlg({'Window size','Overlaping (%)', 'Save full PCD', 'Solver type','Account for artifacts'},'pop_inverseSolution',1,{num2str(round((40/1000)*EEG.srate)),'25', 'bsbl', 'true', 'true'});
    if isempty(answer)
        return;
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
if isempty(solver)
    solver = PEB(H, Delta, blocks);
else
    try
        if sum((solver.H(:) - H(:)).^2) + sum((solver.Delta(:) - Delta(:)).^2) + sum((solver.Blocks(:) - blocks(:)).^2) ~=0
            solver = PEB(H, Delta, blocks); 
        end
    catch ME
        disp(ME)
        solver = PEB(H, Delta, blocks);
    end
end
solver.init;
options = solver.defaultOptions;
options.verbose = false;
if strcmpi(solverType,'loreta')
    options.doPruning = false;
end
EEG.data = double(EEG.data);
Nroi = length(hm.atlas.label);

% Allocate memory
X = allocateMemory([Nx, EEG.pnts, EEG.trials]);
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

windowSize=max([5,windowSize]);
overlapWin = round(windowSize*overlaping);
stepWin = windowSize-overlapWin;
smoothing = hanning(overlapWin*2)';
smoothing = smoothing(1:end/2);

prc_5 = round(linspace(1,EEG.pnts,30));
iterations = 1:stepWin:EEG.pnts-stepWin;
prc_10 = iterations(round(linspace(1,length(iterations),10)));

logE = zeros([length(1:windowSize:EEG.pnts),EEG.trials]);
gamma = zeros([solver.Ng,length(1:windowSize:EEG.pnts),EEG.trials]);
indGamma = EEG.times(1:windowSize:EEG.pnts);

% Perform source estimation
fprintf('PEB source estimation...\n');

for trial=1:EEG.trials
    fprintf('Processing trial %i of %i...',trial, EEG.trials);
    c = 1;
    % for k=1:halfWindow:EEG.pnts
    for k=1:stepWin:EEG.pnts
        loc = k:k+windowSize-1;
        loc(loc>EEG.pnts) = [];
        if isempty(loc), break;end
        if length(loc) < windowSize
            [X(:,loc(1):EEG.pnts,trial),~,~,gamma(:,c,trial), logE(c,trial)] = solver.update(EEG.data(:,loc(1):end,trial), [],[],options);
            break;
        end
        
        
        % Source estimation
        [Xtmp,~,~,gamma(:,c,trial), logE(c,trial)] = solver.update(EEG.data(:,loc,trial),[],[],options);
        
        % Stitch windows
        if k>1 && windowSize > 1
            X(:,loc(1:overlapWin),trial) = bsxfun(@times, Xtmp(:,1:overlapWin), smoothing) + bsxfun(@times,X(:,loc(1:overlapWin),trial), 1-smoothing);
            X(:,loc(overlapWin+1:end),trial) = Xtmp(:,overlapWin+1:end);
        else
            X(:,loc,trial) = Xtmp;
        end
        indGamma(c) = loc(end);
        
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
    
    % Compute average ROI time series
    X_roi(:,:,trial) = computeSourceROI(X, indG, trial, P, isVect);
    
    % Data cleaning
    EEG.data(:,:,trial) = cleanData(H, X, indG, trial);
    
    fprintf('\n');
end
EEG.etc.src.act = X_roi;
EEG.etc.src.roi = hm.atlas.label;
EEG.etc.src.gamma = gamma;
EEG.etc.src.indGamma = indGamma;
EEG.etc.src.H = H;
EEG.etc.src.indG = indG;
EEG.etc.src.indV = indV;
EEG.etc.src.logE = logE;
fprintf('done\n');

if saveFull
    try
        EEG.etc.src.actFull = X;
    catch
        EEG.etc.src.actFull = invSol.LargeTensor([Nx, EEG.pnts, EEG.trials], tempname);
        EEG.etc.src.actFull(:) = X(:);
    end
else
    EEG.etc.src.actFull = [];
end
EEG.history = char(EEG.history,['EEG = pop_inverseSolution(EEG, ' num2str(windowSize) ', ' num2str(overlapWin) ,', ''' solverType ''', ' num2str(saveFull) ', ' num2str(account4artifacts) ');']);
disp('The source estimates were saved in EEG.etc.src');
end


%%
function y = cleanData(H, X, indG, trial)
try
    y = H(:, indG)*X(indG,:, trial);
catch
    n = size(X,2);
    y = zeros(size(H,1),n);
    delta = min([1024 round(n/100)]);
    for k=1:delta:n
        ind = k:k+delta-1;
        ind(ind>n) = [];
        y(:,ind) = H(:, indG)*X(indG,ind, trial);
    end
end
end

%%
function x_roi = computeSourceROI(X, indG, trial, P, isVect)
try
    x = X(indG,:, trial);
    if isVect
        x_roi = sqrt(P*(x.^2));
    else
        x_roi = P*x;
    end
catch
    n = size(X,2);
    x_roi = zeros(size(P,1),n);
    delta = min([1024 round(n/100)]);
    for k=1:delta:n
        ind = k:k+delta-1;
        ind(ind>n) = [];
        if isVect
            x_roi(:,ind) = sqrt(P*(X(indG,ind, trial).^2));
        else
            x_roi = P*X(indG,ind, trial);
        end
    end
end
end

%%
function X = allocateMemory(dim)
try
    X = zeros(dim);
catch ME
    disp(ME.message)
    disp('Using a LargeTensor object...')
    X = LargeTensor(dim);
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
