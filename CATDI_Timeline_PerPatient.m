%% Define Directory and Patient Info

projDir = '/Volumes/bcm-neurosurgery-ecog/ECoG_Data';

%allSubs = {'001','002','003','006','008','010'};
allSubs = {'001','002','006','008','010'};

relativeDays = {'Day01','Day02','Day03','Day04','Day05','Day06','Day07','Day08','Day09','Day10'};

%subStartDates = datetime({'10-Mar-2020';'20-Oct-2020';'13-Apr-2021';'08-Feb-2022';'18-Oct-2022';'09-May-2023'}); %Pulled from patient log b/c not every patient has CAT-DIs on each day of their stay
subStartDates = datetime({'10-Mar-2020';'20-Oct-2020';'08-Feb-2022';'18-Oct-2022';'09-May-2023'}); %Pulled from patient log b/c not every patient has CAT-DIs on each day of their stay

%% Compile Patient Data

allSubScores = table();

for i_sub = 1:numel(allSubs)
    %Current subject info
    currSub = allSubs{i_sub};
    fprintf('%s\n',currSub);
    subFile = [projDir filesep 'DBSTRD' currSub 'Datafile' filesep 'BEHAV' filesep 'scores' currSub '.csv'];
    subTab = readtable(subFile);
    subStartDay = subStartDates(i_sub);
    
    %Parse Date-Time
    DTstrings = cellstr(string(datetime(subTab.Datetime)));
    allDates = arrayfun(@(x) DTstrings{x}(1:strfind(DTstrings{x},' ')-1),1:numel(DTstrings),'UniformOutput',false)';
    allTimes = arrayfun(@(x) DTstrings{x}(strfind(DTstrings{x},' ')+1:end),1:numel(DTstrings),'UniformOutput',false)';
    
    %Add info to table
    subTab.subID = repmat({['DBSTRD' currSub]},size(subTab,1),1);
    subTab.Dates = allDates;
    subTab.Times = allTimes;
    subTab.Hours = hour(allTimes);
    subTab.Minutes = minute(allTimes);
    subTab.RelativeDay = cell(size(subTab,1),1);
    
    %Find Relative Day of Stay
    allSubStayDates = dateshift(subStartDates(i_sub),'start','day',0:9);
    for i_subDate = 1:numel(allDates)
        stayDateIDX = find(allSubStayDates == allDates(i_subDate));
        currRelativeDay = relativeDays{stayDateIDX};
        subTab.RelativeDay{i_subDate} = currRelativeDay;
    end
    
    %Make Full Table
    allSubScores = [allSubScores;subTab];
    
end

%% Make figure for each subject

figure('Name','CAT-DI Assessment Timeline','Position',[0 0 1316 numel(allSubs)*1316/12]);
tt = tiledlayout(numel(allSubs),1);
ylabel(tt,'CAT-DI Score');
xlabel(tt,'Time');

for i_sub = 1:numel(allSubs)
    currSub = allSubs{i_sub};
    currSubLabel = ['DBSTRD' currSub];
    currData = allSubScores(contains(allSubScores.subID,currSubLabel),:);
    
    nexttile;
    title(currSubLabel);
    %xlabel('Time');
    %ylabel('CAT-DI Score');
    
    xlim([0,24*9 + .5]);
    ylim([min(currData.Severity)-5,max(currData.Severity)+5]);
    xticks((([1:1:9]-1)*24)+12);
    xticklabels(relativeDays(1:end-1));
    hold on;
    
    for i_day = 1:numel(relativeDays)
        dayLabel = relativeDays{i_day};
        dayData = currData(contains(currData.RelativeDay,dayLabel),:);
        
        
        xline(i_day*24,'LineWidth',2);
        
        if ~isempty(dayData)
            dayDataCoordinates = (i_day - 1)*24 + dayData.Hours + dayData.Minutes./60;
            scatter(dayDataCoordinates,dayData.Severity,30,'filled','MarkerFaceColor',[58 83 164]./256,'MarkerEdgeColor','k','LineWidth',1.5);
        end
        
        
    end
    
end
set(gcf,'Color','white');
set(gcf,'renderer','Painters');
saveas(gcf,'/Volumes/bcm-neurosurgery-ecog/LABS/Bijanki/Bijanki_Lab/PEOPLE/Maddie/Manuscripts/CATDI_Slope_Biomarker(SUBMITTED_PLOSONE)/CODE/Figures/CATDI_TIMELINE1.eps','epsc');









