%CATDI-Slope biomarker shuffle control check
allSubs = {'DBSTRD001','DBSTRD002','DBSTRD006','DBSTRD008','DBSTRD010'};

projDir = '/Users/mmocchi/Desktop/Carl_Code';
slopesFile = [projDir filesep 'Group_Slope_Fits.mat'];
allSubsSlopes = load(slopesFile);

numPerms = 1000;

for i_sub = 1:numel(allSubs)
    subID = allSubs{i_sub};
    subSlopes = allSubsSlopes.slope_fits{i_sub};
    
    subCatFile = [projDir filesep sprintf('%s_CATDI_Scores.xlsx',subID)];
    subCatTab = readtable(subCatFile);
    subCats = subCatTab.DEP;
    
    goodRunIDX = find(sum(subSlopes,2) ~= 0);
    goodSlopes = subSlopes(goodRunIDX,:);
    goodCats = subCats(goodRunIDX);
    
    pseudoCorrs = NaN(size(subSlopes,2),numPerms);
    obsCorrs = NaN(size(subSlopes,2),2);
    
    tic;
    for i_chan = 1:size(subSlopes,2)
        fprintf('Channel %d\n',i_chan);
        chanSlope = goodSlopes(:,i_chan)';
        
        [obsCorrs(i_chan,1),obsCorrs(i_chan,2)] = corr(chanSlope',goodCats);
        
        [~,randSlopeIDX] = sort(rand(numPerms,size(chanSlope,2)),2);
        
        randChanSlope = arrayfun(@(x) chanSlope(randSlopeIDX(x,:)),[1:numPerms],'UniformOutput',false)';
        randChanSlope = cat(1,randChanSlope{:});
        
        corrDistr = arrayfun(@(x) corr(randChanSlope(x,:)',goodCats),[1:numPerms],'UniformOutput',false)';
        pseudoCorrs(i_chan,:) = cat(1,corrDistr{:});
    end
    toc
    pseudoPerc = prctile(pseudoCorrs,[5 95],2);
    sigChans = obsCorrs(:,1)<=pseudoPerc(:,1) | obsCorrs(:,1) >= pseudoPerc(:,2);
    sigCorrCoeff = obsCorrs(:,2)<=0.05;
    
    falseNeg{i_sub} = ~sigCorrCoeff & sigChans;
    falsePos{i_sub} = sigCorrCoeff & ~sigChans;
end