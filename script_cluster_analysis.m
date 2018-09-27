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
Y = tsne(M','Distance','correlation','NumDimensions',2,'Verbose',1,'Options',struct('MaxIter',1000,'OutputFcn',[],'TolFun',1e-10),'Perplexity',50);

%% K_Means
nc = 15;
[L, C] = kmeans(M',nc,'distance','sqeuclidean', 'Replicates',9);
C = C';

%% Save C L M trainSet and testSet
save(fullfile(results,'cluster.mat'),'C','L','M','trainSet','testSet');

%% Make figure
color = jet(nc);
fig = figure('Position',[243    89   744   700]);
scatter(Y(:,1),Y(:,2),15,color(L,:),'filled');
axis equal
ax_tsne = findobj(fig,'type','axes');
ax_tsne.Position = [0.0681    0.0326    0.8651    0.8994];
hold(ax_tsne, 'on')
grid(ax_tsne, 'on')
set(ax_tsne, 'box','on')
mx = 70;
xlim([-1 1]*mx)
ylim([-1 1]*mx)

ind = reshape(1:nc,3,[]);
figure;
for j=1:nc/3
    for i=1:3
        ax = subplot(3,nc/3,ind(i,j));
        topoplot(C(:,ind(i,j)),refChanlocs,'electrodes','off');
        % topoplot(mean(M(:,L==ind(i,j)),2),refChanlocs,'electrodes','off');
        title(num2str(ind(i,j)))
        colormap(ax,bipolar(256,0.8));
        axis on
        ax.Position(3:4) = [0.1 0.1];
        copyobj(ax,fig);
        mu = [median(Y(L==ind(i,j),1)),median(Y(L==ind(i,j),2))]+1;
        hl = line(ax_tsne,mu(1)*[1 2],mu(2)*[1 2],'Color',color(ind(i,j),:),'LineWidth',2);
    end
end

%% Save figure 
fig.PaperUnits = 'points';
fig.PaperPosition = [0 0 500 500];
print(fig, fullfile(figsFolder,'fig_cluster_labeled.eps'), '-depsc','-r600','-opengl')
