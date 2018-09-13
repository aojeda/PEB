function EEG = pop_inverseSolution(EEG, windowSize, overlaping, solverType, saveFull, postprocCallback)
if nargin < 1, error('Not enough input arguments.');end
if nargin < 5
    answer = inputdlg({'Window size','Overlaping (%)', 'Save full PCD', 'Solver type'},'pop_inverseSolution',1,{num2str((40/1000)*EEG.srate),'50', 'bsbl', 'true'});
    if isempty(answer)
        error('Not enough input arguments.');
    else
        windowSize = str2double(answer{1});
        if isempty(windowSize)
            disp('Invalid input for windowSize parameter, we will use the default value.')
            windowSize= (40/1000)*EEG.srate;
        end
        overlaping = str2double(answer{2});
        if isempty(overlaping)
            disp('Invalid input for overlaping parameter, we will use the default value.')
            overlaping= 50;
        end
        solverType = lower(answer{3});
        if ~any(ismember({'loreta','bsbl'},solverType))
            disp(['Solver ' solverType 'is not available, we will use bsbl.']);
            solverType = 'bsbl';
        end
        saveFull = str2num(lower(answer{4})); %#ok
        if isempty(saveFull)
            disp('Invalid input for saveFull parameter, we will use the default value.')
            saveFull= true;
        end
    end
end
overlaping = overlaping/100;
if nargin < 6, postprocCallback = [];end

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
if exist('Artifact_dictionary.mat','file')
    load('Artifact_dictionary.mat');
    [H, Delta, blocks, indG, indV] = buildAugmentedLeadField(hm, A, chanlocs);
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
        X = invSol.LargeTensor([Nx, EEG.pnts, EEG.trials]);
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

delta = ceil(windowSize*(1-overlaping));
overlap = windowSize-delta;
prc_5 = round(linspace(1,EEG.pnts,30));
iterations = 1:delta:EEG.pnts-windowSize;
prc_10 = iterations(round(linspace(1,length(iterations),10)));

windowSize=max([5,windowSize]);
smoothing = hanning(overlap*2);
smoothing = smoothing(1:overlap)';

gamma = zeros([solver.Ng,length(1:delta:EEG.pnts),EEG.trials]);
time_g = EEG.times(1:delta:EEG.pnts);

% Perform source estimation
fprintf('PEB source estimation...\n');

for trial=1:EEG.trials
    fprintf('Processing trial %i of %i...',trial, EEG.trials);
    c = 1;
    for k=1:delta:EEG.pnts
        loc = k:k+windowSize-1;
        loc(loc>EEG.pnts) = [];
        if isempty(loc), break;end
        if length(loc) < windowSize
            X(:,loc(1):end,trial) = solver.update(EEG.data(:,loc(1):end,trial), [],[],options);
            X_roi(:,loc(1):end,trial) = computeSourceROI(X(indG,loc(1):end,trial),P,isVect);
            break;
        end
        
        % Source estimation
        [Xtmp,~,~,gamma(:,c,trial)] = solver.update(EEG.data(:,loc,trial),[],[],options);
        
        % Stitch windows
        if k>1 && windowSize > 1
            X(:,loc(1:overlap),trial) = bsxfun(@times, Xtmp(:,1:overlap), smoothing) + bsxfun(@times,X(:,loc(1:overlap),trial), 1-smoothing);
            X(:,loc(overlap+1:end),trial) = Xtmp(:,overlap+1:end);
        else
            X(:,loc,trial) = Xtmp;
        end
        
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
    fprintf('\n');
end
fprintf('Done!\n');
EEG.etc.src.act = X_roi;
EEG.etc.src.roi = hm.atlas.label;
EEG.etc.src.gamma = gamma;
EEG.etc.src.time_g = time_g;
EEG.data = EEG.data;
EEG.etc.src.H = H;
EEG.etc.src.indG = indG;
EEG.etc.src.indV = indV;
if saveFull
    try
        EEG.etc.src.actFull = X;
    catch
        EEG.etc.src.actFull = invSol.LargeTensor([Nx, EEG.pnts, EEG.trials], tempname);
    end
else
    EEG.etc.src.actFull = [];
end
EEG.history = char(EEG.history,['EEG = pop_inverseSolution(EEG, ' num2str(windowSize) ', ' num2str(overlaping) ,', ' solverType ', ' num2str(saveFull) ');']);
disp('The source estimates were saved in EEG.etc.src');
end

function x_roi = computeSourceROI(x,P,isVect)
if isVect
    x_roi = sqrt(P*(x.^2));
else
    x_roi = P*x;
end
end