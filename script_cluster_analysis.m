%% Set paths
dataFolders = {'/path/to/study/folder_1', '/path/to/study/folder_2','/path/to/study/folder_n'};
figsFolder = '/path/to/folder/where/figures/will/be/saved';
referenceSetFile = '/path/to/set/file/whose/chanlocs/will/be/used/as/reference/to/coregisted/all/ICs';
results = '/path/to/results/folder';

%% Create train and test file sets
files = [];
for k=1:length(dataFolders)
    tmp = pickfiles(dataFolders{k},'.set');
    files = char(files, tmp);
end
files(1,:) = [];
n = size(files,1);

%% Co-register with Colin27
for subject=1:n
    file = deblank(files(subject,:));
    [filePath,fileName,ex] = fileparts(file);
    if exist(fullfile(filePath,[fileName '_Colin27',ex]),'file')
        continue;
    end
    EEG = pop_loadset(file);
    
    % Remove neck channels
    try         %#ok
        neck = [];
        label = {EEG.chanlocs.labels};
        for ch=1:EEG.nbchan
            if ~isempty(strfind(label{ch},'n'))
                if ~isempty(str2double(label{ch}(2:end)))
                    neck(end+1) = ch;
                end
            end
        end
        EEG = pop_select(EEG,'nochannel',neck);
    end
    
    % Co-register with template
    EEG = pop_forwardModel(EEG);
    
    % Save set
    movefile(EEG.etc.src.hmfile,fullfile(filePath,[fileName '_Colin27.mat']));
    EEG.etc.src.hmfile = fullfile(filePath,[fileName '_Colin27.mat']);
    pop_saveset(EEG, 'filename', [fileName '_Colin27'], 'filepath', filePath, 'savemode', 'onefile'); 
end

%% Make train ans test sets
files = [];
for k=1:length(dataFolders)
    tmp = pickfiles(dataFolders{k},'_Colin27.set');
    files = char(files, tmp);
end
files(1,:) = [];
n = size(files,1);
permutation = randperm(n,n);
files = files(permutation, :);
trainSet = files(randperm(n,round(0.8*n)),:);
testSet = setdiff(files,trainSet,'rows','stable');

%% Co-register ICs
n = size(trainSet,1);
M = [];
template = headModel.loadDefault;
for subject=1:n
    file = deblank(files(subject,:));
    EEG = pop_loadset(file);
    hm = headModel.loadFromFile(EEG.etc.src.hmfile);

    % Interpolate ICs in the channel space of the template
    for ic = 1:size(EEG.icawinv,2)
        F = scatteredInterpolant(hm.channelSpace,EEG.icawinv(:,ic));
        M = [M F(template.channelSpace)];     %#ok
        % Plot IC on the skin of the template to make sure that the interpolation went well 
        % template.plotOnModel(randn(5003,1),M(:,end)); 
    end
end
chanlocs = template.makeChanlocs;


%% T-SNE
Y = tsne(M','Distance','correlation','NumDimensions',2,'Verbose',1,'Options',struct('MaxIter',500,'OutputFcn',[],'TolFun',1e-10),'Perplexity',50);
indEMG = sqrt(sum(Y.^2,2))>40;

%% K_Means
nc = 15;
[L, C] = kmedoids(M',nc,'distance','sqeuclidean', 'Replicates',9);
C = C';

% [L2,C2] = kmeans(M(:,indEMG)',11,'distance','sqeuclidean', 'Replicates',9);
% C2 = C2';


%% Save C L M trainSet and testSet
save(fullfile(results,'cluster.mat'),'C','L','M','trainSet','testSet');

%% Make figure

close all
color = parula(nc);
fig = figure('Position',[243    89   744   700]);
% scatter(Y(:,1),Y(:,2),15,color(L,:),'filled','Marker','o','MarkerFaceAlpha',0.5);
scatter(Y(~indEMG,1),Y(~indEMG,2),15,color(L(~indEMG),:),'filled','Marker','o','MarkerFaceAlpha',0.5);
hold on
colorEMG = [0.5 0.5 0.5];
scatter(Y( indEMG,1)*1,Y(indEMG,2)*1,15,colorEMG,'filled','Marker','o','MarkerFaceAlpha',0.5);
axis equal
ax_tsne = findobj(fig,'type','axes');
ax_tsne.Position = [0.0681    0.0326    0.8651    0.8994];
hold(ax_tsne, 'on')
grid(ax_tsne, 'on')
set(ax_tsne, 'box','on')
mx = 60;
xlim([-1 1]*mx)
ylim([-1 1]*mx)

figure;
M1 = M(:,~indEMG);
ang = linspace(0,360*(nc-1)/nc,nc)*pi/180;%+5*pi/6;
for i=1:nc
    [ax_x,ax_y] = pol2cart(ang(i),0.65);
    ax = axes('Position',[ax_x/2+0.5-0.05,ax_y/2+0.5-0.05,0.1,0.1]);
    %topoplot(C(:,i),chanlocs,'electrodes','off');
    topoplot(mean(M1(:,L(~indEMG)==i),2),chanlocs,'electrodes','off');
    title(num2str(i))
    colormap(ax,bipolar(256,0.8));
    axis(ax,'on');
    ax.Position(3:4) = [0.1 0.1];
    ax = copyobj(ax,fig);
    axis(ax,'off');
end


% fig = figure('Position',[243    89   744   700]);
% scatter(Y2(:,1),Y2(:,2),15,colorEMG(L2,:),'filled','MarkerFaceAlpha',0.75);
% figure

M2 = M(:,indEMG);
nc2 = size(C2,2);
ang = linspace(0,360*(nc2-1)/nc2,nc2)*pi/180+pi/2;
for i=1:nc2
    [ax_x,ax_y] = pol2cart(ang(i),0.85);
    ax = axes('Position',[ax_x/2+0.5-0.05,ax_y/2+0.5-0.05,0.1,0.1]);
    %topoplot(C2(:,i),chanlocs,'electrodes','on');
    topoplot(median(M2(:,L(indEMG)==i),2),chanlocs,'electrodes','off');
    title(num2str(i))
    colormap(ax,bipolar(256,0.8));
    axis(ax,'on');
    ax.Position(3:4) = [0.1 0.1];
    ax = copyobj(ax,fig);
    axis(ax,'off');
end

%% Save figure
fig.PaperUnits = 'points';
fig.PaperPosition = [0 0 500 500];
print(fig, fullfile(figsFolder,'fig_cluster_labeled.eps'), '-depsc','-r600','-opengl')
