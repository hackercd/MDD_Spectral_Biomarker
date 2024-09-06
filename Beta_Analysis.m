%% Set up Paths
addpath('H:\Process')

% INITIALIZE
clear bn
for pi=1:5
    bn{pi}=[];
    patids = {'DBSTRD001','DBSTRD002','DBSTRD006','DBSTRD008','DBSTRD010'};
    patid=patids{pi};
    datadir=['H:\DBSTRD\SEEG\' patid '\Data\'];
    root='H:\DBSTRD\SEEG\';
    Subj_Data_root = sprintf('%s/%s/Data/',root,patid);
    switch pi
        case {1,3,4,5}
            CATfiles = dir(['H:\DBSTRD\SEEG\' patid '\Data\*CAT*\*ns5']);
        case 2
            CATfiles = dir(['H:\DBSTRD\SEEG\' patid '\Data\*CAT*ns3']);
    end
    %     bn{pi}= erase({CATfiles.name},".ns3");
    CATraw{pi} = readtable(['H:\DBSTRD\SEEG\' patid '_CATDI_Scores.xlsx']);
    CAT = CATraw{pi}.DEP;
    CATstore{pi}=CAT;
    for i = 1:length(CATraw{pi}.SEEGBlock)
        part=strsplit(CATraw{pi}.SEEGBlock{i},' ');
        bn{pi}{i}=part{end};
    end
end
% Load metadata
load('Group_CAT_Spectra_RefMedian.mat')

% Set up electrode grouping (multi-patient)
for pi = 1:5
    key_cat{pi}=[];
    for i = find(strcmp(electrodes_all{pi}.Type,'sEEG'))'
        x=split(electrodes_all{pi}.Label(i),'-');
        key_cat{pi} =  [key_cat{pi};string(x{2}(1:end-2))];
    end
    unique(key_cat{pi})
end

keys{1} = {'-ACC','-Amy','-OF','-mOF','-vmPF'};
keys{2} = {'-ACC','-Amy','-OF','-mOF','-mPF'};
keys{3} = {'-ACC','-Amy','-OF','-mOF','-vmP'};
keys{4} = {'-ACC','-Amy','-OF','-mOF','-mPF'};
keys{5} = {'-ACC','-Amy','-OFC','-mOF','-mPF'};
master_keys = {'LACC','LAmy','LOFC','LmOF','LmPF','RACC','RAmy','ROFC','RmOF','RmPF','LSCC','RSCC','LVCVS','RVCVS'};
master_regions = {'ACC','Amy','OFC','mOF','mPF'};

% Set up table
varNames = [{'Index'},{'Electrodes'}, patids];
emptyVars = cell(1, numel(varNames));
electrodes_coreg = table(emptyVars{:}, 'VariableNames', varNames);

k=0;
for i = 1:10 %just SEEG
    for j = 1:16
        k=k+1;
        electrodes_coreg.Electrodes{k}=[master_keys{i} sprintf('%02d',j)];
    end
end

for pi = 1:5
    key_cat{pi}=[];
    for i = find(strcmp(electrodes_all{pi}.Type,'sEEG'))'
        index = cellfun(@(key) contains(electrodes_all{pi}.Label{i}, key), keys{pi});
        leadstr = master_regions{index};
        side = electrodes_all{pi}.Hemisphere{i}(1);
        num=(electrodes_all{pi}.Label{i}(end-1:end));
        matchstr = [side leadstr num];
        foundind = contains(electrodes_coreg.Electrodes,matchstr);
        if sum(foundind)~=1
            error(sprintf('found %d matches\n',sum(foundind)))
        end
        electrodes_coreg.(patids{pi})(foundind)=i;
    end
end
electrodes_coreg.Index=[1:size(electrodes_coreg,1)]';

% Load slope fits / FOOOFed beta power
load('Group_FOOOF_Results_15_45.mat')

% Define valid blocks/electrode
for pi=1:5
    patid=patids{pi}
    block_skip{1} = [36];
    block_skip{2} = [2,19];
    block_skip{3} = [30,38];
    block_skip{4} = [44];
    block_skip{5} = [24 66,67];
    switch pi
        case 1
            CARignore=[1:16,23,24,53,54,70:71,88:103];
            CARignore_all{pi}=CARignore;
            CARignore=unique([CARignore+1 CARignore]);
        case 2
            CARignore=[66,67,1:8,22:29,72:79,94:101];
            CARignore=[CARignore 69 144 55 11 137];
            CARignore_all{pi}=CARignore;
        case 3
            CARignore=[1:32 33 48 64 80 107 108 122 138 152 153 154 169 170 184];
            CARignore_all{pi}=CARignore;
        case 4
            CARignore=[1:32 72 79 104 111 113 118 132 161 162 163 164 180];
            CARignore_all{pi}=CARignore;
        case 5 %this will be in ns3 ordering until .csv
            CARignore=[1:32];
            CARignore_all{pi}=CARignore;
    end
end

%%
fooof_fits_neg = cellfun(@(x) -x, fooof_fits, 'UniformOutput', false);
SpecStats = {fooof_fits_neg,beta_lows_orig,beta_lows,beta_his_orig,beta_his};
clear r_FvC_VM
for comparator = 1:length(SpecStats)
    clear r_FvCs p_FvCs
    VM = [33:38, 49:52, 65:69, 113:118, 129:132, 145:148];
    ROIs(VM) =   1;
    r_FvC_VM{comparator} = [];
    
    for pi=1:5;
        valid_blocks{pi}=1:size(fooof_fits{pi},1);
        valid_blocks{pi}(block_skip{pi}(block_skip{pi}<size(fooof_fits{pi},1)))=[];
        c1=CATraw{pi}.DEP(valid_blocks{pi});

        ve=1:size(electrodes_all{pi},1);
        ve([CARignore_all{pi}])=[];
        for j = 1:size(fooof_fits{pi},2)
            f1 = SpecStats{comparator}{pi}(valid_blocks{pi},j);
            [r_FvCs{pi}(j),p_FvCs{pi}(j)] = corr(c1,f1);
        end
        
        
        ROIp = nonzeros(electrodes_coreg{ROIs==1,patids{pi}});
        ROIp=ROIp(~ismember(ROIp,CARignore_all{pi}));
        f1 = median(SpecStats{comparator}{pi}(valid_blocks{pi},ROIp),2);
        [r_FvC_means{comparator,1}(pi),p_FvC_means{comparator,1}(pi)] = corr(c1,f1);

        %r_FvC_means{comparator,1}(pi)=median(r_FvCs{pi}(ROIp)');
    end

    
    for pi =1:5
        ROIp = nonzeros(electrodes_coreg{ROIs==1,patids{pi}});
        ROIp=ROIp(~ismember(ROIp,CARignore_all{pi}));
        r_FvC_VM{comparator} = [r_FvC_VM{comparator}; r_FvCs{pi}(ROIp)'];
    end
end

r_FvC_ROI_cat = vertcat(r_FvC_VM{:});
grouping = repelem(1:numel(r_FvC_VM), cellfun('length', r_FvC_VM));


figure(1),clf
violinplot(r_FvC_ROI_cat,grouping,'ViolinColor',.5*ones(5,3),'BoxColor',[0 0 0],'linewidth',0.3)
set(gca,'xticklabels',{'Aperiodic Slope','Low Beta (Total Power)','Low Beta (Osc. only)','High Beta (Total Power)','High Beta (Osc. only)'},'fontsize',16)
ylabel('Depression Severity Correlation','FontSize',18),ylim([-0.85,.85])
line(xlim,[0 0],'color','k')
[h,p,ci,stats] = ttest2(r_FvC_ROI_cat(grouping==1),r_FvC_ROI_cat(grouping==2))
for i = 1:length(SpecStats)
    x=rand(5,1);
    for pi = 1:5
        plot(i+x(pi)*.4-.2,r_FvC_means{i}(pi),'o','markersize',15,'markerfacecolor','w','markeredgecolor','k')
        text(i+x(pi)*.4-.2,r_FvC_means{i}(pi),num2str(pi),'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
end

%%
for i = 1:length(SpecStats)
    [h,p,ci,stats] = ttest(r_FvC_means{i});
%     [p,h,stats] = signrank(r_FvC_means{i});
    p
    
end










%% Set up group matrices
clear group_coords
for pi = 1:5
    for i = 1:size(electrodes_coreg,1)
        if ~~electrodes_coreg.(patids{pi})(i)
            group_r_F(pi,i) = r_FvCs{pi}(electrodes_coreg.(patids{pi})(i));
            group_coords(pi,i,1) = electrodes_all{pi}.Coord_x(electrodes_coreg.(patids{pi})(i));
            group_coords(pi,i,2) = electrodes_all{pi}.Coord_y(electrodes_coreg.(patids{pi})(i));
            group_coords(pi,i,3) = electrodes_all{pi}.Coord_z(electrodes_coreg.(patids{pi})(i));
        else
            group_r_F(pi,i) = nan;group_coords(pi,i,1:3)=nan;
        end
    end
end
all_coords = group_coords;

mean_coords=squeeze(mean(group_coords,1));
std_coords=squeeze(mean(group_coords,1));
missing=~~sum(isnan(group_coords(:,:,1)),1);
for i=1:10
    clear interpCoords
    indices=(1:16)+(i-1)*16;
    coords= mean_coords(indices,:);
    x=find(~missing(indices));
    y=coords(~missing(indices),:);
    x1=find(missing(indices));
    for d=1:3
        interpCoords(:,d) = interp1(x, y(:,d), x1, 'linear', 'extrap');
    end
    mean_coords(x1+(i-1)*16,:) = interpCoords;
end

%% Exemplar spectra
pi=1;roii=1;


figure(2),clf
for pi = 1:5

end
%%
clear x fit_x xlog
pi=3
eoi=77 %pi 1 115  % pi 1 146 - high knee
[~,fxx] = pwelch(rand(10000,1),2000,[],[],2000);
f_log = logspace(log10(.1), log10(200), 200)';

f_fit_bins = [20 45];
foi=f_log>f_fit_bins(fi,1)&f_log<f_fit_bins(fi,2);

for i= valid_blocks{pi}
    x(:,i) = signals_pxx{pi}{i}(:,eoi);
end
for e = 1:size(x,2)
    xlog(:,e) = 10.^interp1(fxx,log10(x(:,e)),f_log);
    fit_x(:,e) = polyfit(log10(f_log(foi)),log10(xlog(foi,e)), 1);
end
C=CATstore{pi}(valid_blocks{pi});
ftx=fit_x(1,valid_blocks{pi});
figure(1),clf
subplot 211
plot(ftx,C,'o')
subplot 212,hold on
[~,mini]=min(fit_x(1,:));
[~,maxi]=max(fit_x(1,:));

plot(fxx,signals_pxx{pi}{mini}(:,eoi),'k--')
plot(fxx,signals_pxx{pi}{maxi}(:,eoi),'k')
set(gca,'yscale','log','xscale','log')
%         axis([3 130 1e-2 5e2])
xlim([3 120])
set(gca,'xtick',[1 10 20 45 100])
legend({'More Depressed','Less Depressed'})

%% Plot 3d text for pi, field
% field = 'ROI_vis';
field = 'SortedIndex';
% pi=3;
doplotsurf=1;
if 1, [vl,fl]=read_surf('S:\Imaging Data\Surfaces\lh.pial');
    [vr,fr]=read_surf('S:\Imaging Data\Surfaces\rh.pial');end
figure(12),set(gcf,'color','w'),clf
Y=r_FvCs;
cm=colormap(GenerateCMap('RedGreyBlue',.15));
cl=[-.8 .8];
for k=1%1:length(Y)
    clf,hold on,axis equal off vis3d, view([-180,0])
    patch('faces',fl+1,'vertices',vl,'facealpha',.1,'facecolor',[0.7 .7 .7],'edgealpha',.00,'FaceLighting','phong','EdgeLighting','phong')%,shading interp
    patch('faces',fr+1,'vertices',vr,'facealpha',.1,'facecolor',[0.7 .7 .7],'edgealpha',.00,'FaceLighting','phong','EdgeLighting','phong')
    light
    c_int = floor((Y{k}-cl(1))/(cl(2)-cl(1))*255)+1;
    c = squeeze(ind2rgb(c_int,cm));
    for i = 1:size(electrodes_all{k},1)
        x=electrodes_all{k}.MNI305_x(i);
        y=electrodes_all{k}.MNI305_y(i);
        z=electrodes_all{k}.MNI305_z(i);
%         text(x,y,z,electrodes_all{k}.(field){i});
        text(x,y,z,num2str(electrodes_all{k}.(field)(i)));
        plot3(x,y,z,'o','markerfacecolor',c(i,:),'markeredgecolor','k','markersize',12);
        %         if ~mod(i,16)
        
        %         end
    end
end
