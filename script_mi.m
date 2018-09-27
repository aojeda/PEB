clear all
dataFolder = '/data/ESI_ICA';
load('/data/ESI_ICA/cluster.mat','testSet');
%%

ntrials = [1 2 5 10 20 40 80];
Ntrials = length(ntrials);
MI = zeros(size(testSet,1),4, Ntrials);

hwait = waitbar(0);
for k=3:size(testSet,1)
    file = deblank(testSet(k,:));
    if ~exist(file,'file')
        [~,fname,ex] = fileparts(file);
        file = fullfile(dataFolder,[fname, ex]);
    end
    EEG = pop_loadset(file);
    EEG = pop_forwardModel(EEG, headModel.getDefaultTemplateFilename, [0.33, 0.022, 0.33], 1,1);
    EEG = pop_eegfiltnew(EEG, [],1,1650,1,[],0);
    EEG = pop_resample(EEG, 250);
    EEG = pop_reref( EEG, []);
    for t=1:Ntrials
        EEG_epoch = pop_epoch(EEG, unique({EEG.event.type}), [-2 2]);
        EEG_epoch = pop_select(EEG_epoch, 'trial', randperm(EEG_epoch.trials,min([ntrials(t) EEG_epoch.trials])));
        
        EEG_pebp = pop_inverseSolution(EEG_epoch, 20, 25,'bsbl',false, true);
        EEG_pebc = pop_inverseSolution(EEG_epoch, 20, 25,'bsbl',false, false);
        EEG_epoch = pop_runica(EEG_epoch, 'extended',1);    
        
        MI1 = minfo(EEG_epoch.data);
        MI1 = MI1 - diag(diag(MI1));
        MI(k,1,t) = sum(sum(MI1))/numel(MI1);
    
        MI2 = minfo(EEG_epoch.icaact);
        MI2 = MI2 - diag(diag(MI2));
        MI(k,2,t) = sum(sum(MI2))/numel(MI2);
        
        MI3 = minfo(EEG_pebc.etc.src.act);
        MI3 = MI3 - diag(diag(MI3));
        MI(k,3,t) = sum(sum(MI3))/numel(MI3);
        
        MI4 = minfo(EEG_pebp.etc.src.act);
        MI4 = MI4 - diag(diag(MI4));
        MI(k,4,t) = sum(sum(MI4))/numel(MI4);
    end
    clear EEG_epoch EEG_pebp EEG_pebc;
    waitbar(k/size(testSet,1), hwait);
end
close(hwait);
%%
fig = figure;
MIR = bsxfun(@minus,MI(:,1,:),MI(:,[2 4],:));
ax = subplot(121);
bar(categorical({'4' '8' '20' '40' '80' '160' '320'},{'4' '8' '20' '40' '80' '160' '320'}),squeeze(median(MIR,1))');
xlabel('Data segment legth (sec)');
ylabel('Mutual information reduction')
legend({'Infomax ICA','PEB+'})
ax.XAxis.FontSize = 7;
ax.YAxis.FontSize = 7;
ax.YLim = [0.12 0.1825]; % tinv(0.95,27-1);

ax = subplot(122);
MIR = squeeze(MI(:,2,:) - MI(:,4,:));
boxplot(MIR,'Labels',{'4' '8' '20' '40' '80' '160' '320'});grid on
xlabel('Data segment legth (sec)');
ylabel('Mutual information reduction (PEB+ vs Infomax ICA)')
ax.YLim = [-.05 .05];
ax.XAxis.FontSize = 7;
ax.YAxis.FontSize = 7;
%%
figsFolder = '/home/ale/Dropbox/ESI_ICA/';
fig.PaperUnits = 'points';
fig.PaperPosition = [0 0 500 200];
print(fig, fullfile(figsFolder,'fig_mir.eps'), '-depsc','-r600','-painters')
%%
fig = figure;
plot([MI_eeg, MI_ica MI_pebc, MI_pebp])