%% INITIALIZE
% Acquire data file strucutre corresponding to CAT test scores
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

%% Acquire CAT times for crossreference to CAT score sheet
% Only run once per patietn
% Strings transferred to spreadsheet for manual verification
for i = 1:length(CATfiles)
    ns5_root = [CATfiles(i).folder filesep CATfiles(i).name];
    NS3=openNSx(ns5_root,'noread');
    fprintf('%s %s\n',NS3.MetaTags.DateTime,erase(CATfiles(i).name,'.ns5'))
end

%% Convert data, run once per patient
for pi = 1:5
    for i = 1:length(bn{pi})
        clear x
        % Load raw data info for each block
        block_name = bn{pi}{i};
        ns3_root = [CATfiles(i).folder filesep CATfiles(i).name(1:end-3) 'ns3'];
        montageInfo = blackrockMontageInfo_modified(ns3_root,[]);
        
        % Read raw data (read NS5 if NS3 not available)
        if exist(ns3_root)
            NS3=openNSx(ns3_root,'read','p:double');
        else
            %no NS3 - check NS5
            ns3_root = [CATfiles(i).folder filesep CATfiles(i).name(1:end-3) 'ns5'];
            if ~exist(ns3_root)
                error('data file not found');
            end
            NS3=openNSx(ns3_root,'read','p:double');
            disp('decimating 30k -> 2k')
            for j = 1:size(NS3.Data,1)
                x(j,:)=decimate(NS3.Data(j,:),NS3.MetaTags.SamplingFreq/2000);
            end
            NS3.Data=x;
        end
        
        % Define master parameters        
        master_vars.ecog_srate=2000;
        master_vars.compress=2;
        master_vars.badchan = [];
        
        % Give consistent naming for electrode montages
        fprintf('%s %s\n',NS3.MetaTags.DateTime,block_name)
        switch pi
            case {1,2}
                master_vars.originalelecs= montageInfo.sEEG.Contacts.Indices;
            case {3,4,5}
                master_vars.originalelecs= [montageInfo.DBS.Contacts.Indices;montageInfo.sEEG.Contacts.Indices];
        end
        
        % Define montage parameters        
        master_vars.probes=montageInfo.sEEG.Probes.Indices;
        master_vars.nchan = size(master_vars.originalelecs,1);
        master_vars.elecs= [1:master_vars.nchan];
                
        % Create index_csv_to_NS3 - vector indexing electrode metadata to raw data
        electrodes_orig = readtable(sprintf('%s/%s-electrodes.csv',root,patid));
        csvlabels = electrodes_orig.Label;
        NS3labels = {NS3.ElectrodesInfo(master_vars.originalelecs).Label}';
        index_csv_to_NS3 = matchstrings_csvtoNS3(csvlabels,NS3labels,pi);
        master_vars.index_csv_to_NS3=index_csv_to_NS3;
        electrodes = [table([1:size(master_vars.originalelecs,1)]','VariableNames',{'SortedIndex'}),electrodes_orig(master_vars.index_csv_to_NS3,:)];
        
        % Save electrode metadata to master parameters
        master_vars.electrodes=electrodes;
        
        % Save master parameters
        fn=sprintf('%smaster_%s.mat',Subj_Data_root,block_name);
        save(fn,'master_vars');
        
        % Save raw data
        RawVoltage = NS3.Data(master_vars.originalelecs,:);
        RawVoltage = RawVoltage*0.25;
        fn=sprintf('%s\\%s\\Data\\RawVoltage_%s.mat',root,patid,block_name);
        save(fn,'RawVoltage','-v7.3','-nocompression');
    end
end

%% QA electrode label/sorting
% Optional - used to verify metadata matches correct electrode indices
electrodes_orig = readtable(sprintf('%s/%s-electrodes.csv',root,patid));
csvlabels = electrodes_orig.Label;
NS3labels = {NS3.ElectrodesInfo(master_vars.originalelecs).Label}';
master_vars.index_csv_to_NS3=matchstrings_csvtoNS3(csvlabels,NS3labels,pi);
% Check that these line up ---> [NS3labels,csvlabels(master_vars.index_csv_to_NS3)]
electrodes = [table([1:size(master_vars.originalelecs,1)]','VariableNames',{'SortedIndex'}),electrodes_orig(master_vars.index_csv_to_NS3,:)]

%% Compute Spectra 
clear signals_filt signals_filt_mean signals_pxlogx fit_pxx mean_pxx
clear signals_pxx
figure(1)
% Looping over patients
for pi=1:5
    patid=patids{pi}
    % Bad blocks - epochs overly contaminated by noise, patient motion
    block_skip{1} = [36];
    block_skip{2} = [2,19];
    block_skip{3} = [30,38];
    block_skip{4} = [44];
    block_skip{5} = [24 66,67];
    % Electrodes removed from analysis - includes all DBS contacts and SEEG
    % electrodes contaminated by exces
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
    
    valid_blocks{pi}=1:length(bn{pi});
    valid_blocks{pi}(block_skip{pi})=[];
    
    filtfs=[5,7.5;8,12;70,110];
    % continuous estimator bins
    fcs=logspace(log10(20),log10(45),5);
    fhw=5;
    % log spaced f params
    
    for i = 1:length(bn{pi})
        if ismember(i,block_skip{pi}),continue,end
        block_name = bn{pi}{i}
        Subj_Data_root = sprintf('%s/%s/Data/',root,patid);
        fn=sprintf('%smaster_%s.mat',Subj_Data_root,block_name);    load(fn);
        fn=sprintf('%s\\%s\\Data\\RawVoltage_%s.mat',root,patid,block_name);    load(fn);
        signals=RawVoltage';clear RawVoltage;
        fs=master_vars.ecog_srate;
        switch pi
            case 1
                electrodes_orig = readtable(sprintf('%s/%s-electrodes.csv',root,patid));
                electrodes = [table([1:size(master_vars.originalelecs,1)]','VariableNames',{'SortedIndex'}),electrodes_orig(master_vars.originalelecs,:)];
            case 5
                %             electrodes = table([1:size(master_vars.originalelecs,1)]');
                
                electrodes_orig = readtable(sprintf('%s/%s-electrodes.csv',root,patid));
                electrodes = [table([1:size(master_vars.originalelecs,1)]','VariableNames',{'SortedIndex'}),electrodes_orig(master_vars.index_csv_to_NS3,:)];
            otherwise
                electrodes_orig = readtable(sprintf('%s/%s-electrodes.csv',root,patid));
                electrodes = [table([1:size(master_vars.originalelecs,1)]','VariableNames',{'SortedIndex'}),electrodes_orig(master_vars.index_csv_to_NS3,:)];
        end
        
        % Rereference to common median
        allchannels=1:size(electrodes,1);
        RR_channels = allchannels(~ismember(allchannels,CARignore));
        signals_rr = signals-repmat(median(signals(:,RR_channels),2),1,size(signals,2));

        % compute spectra
        [signals_pxx{pi}{i},fxx] = pwelch(signals_rr,fs,[],[],fs);
        electrodes_all{pi}=electrodes;
        signal_lengths{pi}(i)=size(signals_rr,1);
    end
end
% save('Group_CAT_Spectra_Ref.mat','signals_pxx','electrodes_all','block_skip')

%% Electrode region labeling keys - Coregister electrode identities across subjects
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
for i = 1:10 %just SEEG indices
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

%% Set up group coordinate matrices
clear group_coords
for pi = 1:5
    for i = 1:size(electrodes_coreg,1)
        if ~~electrodes_coreg.(patids{pi})(i)
            group_r_F(pi,i) = r_FvCs{pi}(electrodes_coreg.(patids{pi})(i));
            group_coords(pi,i,1) = electrodes_all{pi}.Coord_x(electrodes_coreg.(patids{pi})(i));
            group_coords(pi,i,2) = electrodes_all{pi}.Coord_y(electrodes_coreg.(patids{pi})(i));
            group_coords(pi,i,3) = electrodes_all{pi}.Coord_z(electrodes_coreg.(patids{pi})(i));
        else
            group_r_F(pi,i) = nan;
            group_coords(pi,i,1:3)=nan;
        end
    end
end
all_coords = group_coords;
mean_coords=squeeze(mean(group_coords,1));

%% Plot Group effect
field = 'ROI_vis';pi=1;doplotsurf=1;
% Load surfaces
if 1, [vl,fl]=read_surf('S:\Imaging Data\Surfaces\lh.pial');
    [vr,fr]=read_surf('S:\Imaging Data\Surfaces\rh.pial');end
figure(11),set(gcf,'color','w'),clf
clear Y
Y{1}=mean(group_r_F,'omitnan');
Y{2}=mean(group_r_F,'omitnan')./std(group_r_F,'omitnan');
cls={[-.8 .8],[-5 5]};
cm=colormap(GenerateCMap('RedGreyBlue',.15,[.9 .9 .9]));
skip = [12:16,30:32,110:112,93:96,48,64,80,144,128,158:160]; % ignore group indices out of shaft range
for k=1%1:length(Y)
    cl=cls{k};
    clf,hold on
    h1=patch('faces',fl+1,'vertices',vl,'facealpha',.1,'facecolor',[0.7 .7 .7],'edgealpha',.00,'FaceLighting','phong','EdgeLighting','phong')%,shading interp
    h2=patch('faces',fr+1,'vertices',vr,'facealpha',.1,'facecolor',[0.7 .7 .7],'edgealpha',.00,'FaceLighting','phong','EdgeLighting','phong')
    light
    c_int = floor((Y{k}-cl(1))/(cl(2)-cl(1))*255)+1;
    c = squeeze(ind2rgb(c_int,cm));
    for i = size(electrodes_coreg,1):-1:1
        if ismember(i,skip),continue,end
        if sum(isnan(group_r_F(:,i)))<2
            x=mean_coords(i,1);
            y=mean_coords(i,2);
            z=mean_coords(i,3);
            plot3(x,y,z,'o','markerfacecolor',c(i,:),'markeredgecolor','k','markersize',abs(Y{2}(i))*3+5); % ensure scaling consistent with legend
        end
    end
    cb=colorbar('FontSize',11),caxis(cl);
    axis equal off vis3d, view([-180,0]);
end

%% ROI Violin Plots
clear ROIs
ACC = [1:4, 81:84];
Amy = [17:21, 97:100];
Temp = [24:29, 104:109];
VM = [33:38, 49:52, 65:69, 113:118, 129:132, 145:148];
DL = [7:11, 42:47, 58:63, 73:79, 87:92, 122:127, 137:143,153:157];
ROIs(VM) =   1; ROIs(Amy1) =  2; ROIs(Amy2) =  3; ROIs(ACC) =  4; ROIs(Temp) = 5; for i = 1:5
    r_FvC_ROI{i} = [];
    for pi =1:5
        ROIp = nonzeros(electrodes_coreg{ROIs==i,patids{pi}});
        ROIp=ROIp(~ismember(ROIp,CARignore_all{pi}));
        r_FvC_ROI{i} = [r_FvC_ROI{i}; r_FvCs{pi}(ROIp)'];
    end
end

r_FvC_ROI_cat = vertcat(r_FvC_ROI{:});
grouping = repelem(1:numel(r_FvC_ROI), cellfun('length', r_FvC_ROI));
cs=[1 0 1;0 1 1;1 0 0;0 1 0;0 0 1];
        cs(1,:) = [1 0 1];
%         cs(2,:) = [255,127,0]/255;%[0 1 1];
%         cs(3,:) = [228,26,28]/255;%[1 0 0];
        cs(4,:) = [77,175,74]/255;[0 1 0];
        cs(5,:) = [55,126,184]/255;[0 0 1];
        cs(7,:) = [.9 .9 .9];

        
        cs(3,:) = [255,127,0]/255;%[0 1 1];
        cs(2,:) = [228,26,28]/255;%[1 0 0];
        
figure(1),clf
violinplot(r_FvC_ROI_cat,grouping,'ViolinColor',cs,'BoxColor',[0 0 0],'linewidth',0.3)
set(gca,'xticklabels',{'Ventral PFC','Amygdala','Dorsal PFC','ACC','Temporal'},'fontsize',16)
ylabel('Slope Correlation','FontSize',18),ylim([-0.8,.8])
line(xlim,[0 0],'color','k')
[h,p,ci,stats] = ttest2(r_FvC_ROI_cat(grouping==1),r_FvC_ROI_cat(grouping==2))