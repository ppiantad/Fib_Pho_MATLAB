%% Requires MATLAB 2019 or above for optimal execution
% Requires the following additional files: 
% tblvertcat (compatible w/ R2019b and R2020a): https://www.mathworks.com/matlabcentral/fileexchange/79912-tblvertcat
% TDT MATLAB SDK: https://www.tdt.com/support/matlab-sdk/
% using Scientific Colour-Maps 6.0 (http://www.fabiocrameri.ch/colourmaps.php)
% Generate maximally perceptually-distinct colors: https://www.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors


%% set root directory (where files are located on PC)
root = upper('h:');

which_region = input('Which dataset do you want to analyze?\n 1) Code Test \n 2) BLA-NAcSh \n 3) vmPFC-NAcSh \n 4) D1 \n 5) D2 \n 6) BLA-NAcSh GFP \n 7) vmPFC-NAcSh GFP \n 8) GFP (to NAcSh) ALL \n 9) D1-iOP \n 10) D2-iOP \n 15) D1-eOP \n 16) D2-eOP \n 12) vHPC-NAcSh \n 13) BLA-PL \n','s');

which_sessions = input('Which sessions do you want to analyze from these mice? \n 1) Early RM \n 2) Late RM \n 3) RDT \n 4) Shock Test \n','s');


[animalNames, blockpaths, behavFiles, boris_files,  SLEAP_files, SLEAP_time_range_adjustment, session, implant_side, large_rew_side, whichStreams, whichTTL] = read_mouse_data_v2(root, which_region, which_sessions);


%%
numAnimals=numel(animalNames);
dataStruct_names=animalNames;
timeShift=zeros(1,numAnimals);
Channel_405_name=animalNames;
Channel_465_name=animalNames;
TTL_name=animalNames;
numTrials=animalNames;

ID_Table = cell2table(animalNames);

%make table with all animal info

Table=table(animalNames,blockpaths',behavFiles',boris_files',SLEAP_files', SLEAP_time_range_adjustment', dataStruct_names,timeShift',Channel_405_name,Channel_465_name,TTL_name,whichStreams,whichTTL,numTrials,session', implant_side', large_rew_side');
Headers={'animalNames','blockpath','behavFiles','boris_files','SLEAP_files','SLEAP_time_range_adjustment','dataStruct_names','timeShift','Channel_405_name','Channel_465_name','TTL_name','whichStreams','whichTTL','numTrials','session','implant_side','large_rew_side'};
Table.Properties.VariableNames([1:17])=Headers;


%% Select which session to look at
sessionNum = input('Please enter the # of sessions to analyze per mouse (1, 2, 3, all): \n','s');
% sessionNum = split(sessionNum,",");


% Table = Table(any(Table.session == [str2double(sessionNum(1)],2),:);
% numAnimals = numel(Table.animalNames(:,1));
% ID_Table = cell2table(Table.animalNames(:,1));


if sessionNum == '1'
    idx = Table.session == 1;
    Table = Table(idx,:);
    numAnimals = numel(Table.animalNames(:,1));
    ID_Table = cell2table(Table.animalNames(:,1));
elseif sessionNum == '2'
    idx = Table.session == 2;
    Table = Table(idx,:);
    numAnimals = numel(Table.animalNames(:,1));
    ID_Table = cell2table(Table.animalNames(:,1));
elseif sessionNum == '3'
    idx = Table.session == 3;
    Table = Table(idx,:);
    numAnimals = numel(Table.animalNames(:,1));
    ID_Table = cell2table(Table.animalNames(:,1));
elseif sessionNum == '1,2'
    idx = Table.session == 1 | Table.session == 2;
    Table = Table(idx,:);
    numAnimals = numel(Table.animalNames(:,1));
    ID_Table = cell2table(Table.animalNames(:,1));
elseif sessionNum == '2,3'
    idx = Table.session == 2 | Table.session == 3;
    Table = Table(idx,:);
    numAnimals = numel(Table.animalNames(:,1));
    ID_Table = cell2table(Table.animalNames(:,1));
elseif sessionNum == 'all'
end

%%
impSide = input('Please enter implant side(s) to include (1 = left,  2 = right, all): \n','s');

if impSide == '1'
    idx = Table.implant_side == 1;
    Table = Table(idx,:);
    numAnimals = numel(Table.animalNames(:,1));
    ID_Table = cell2table(Table.animalNames(:,1));
elseif impSide == '2'
    idx = Table.implant_side == 2;
    Table = Table(idx,:);
    numAnimals = numel(Table.animalNames(:,1));
    ID_Table = cell2table(Table.animalNames(:,1));
elseif impSide == 'all'
end

%%
LargeRewSide = input('Please enter Large Rew side to include (1 = left,  2 = right, all): \n','s');

if LargeRewSide == '1'
    idx = Table.large_rew_side == 1;
    Table = Table(idx,:);
    numAnimals = numel(Table.animalNames(:,1));
    ID_Table = cell2table(Table.animalNames(:,1));
elseif LargeRewSide == '2'
    idx = Table.large_rew_side == 2;
    Table = Table(idx,:);
    numAnimals = numel(Table.animalNames(:,1));
    ID_Table = cell2table(Table.animalNames(:,1));
elseif LargeRewSide == 'all'
end



%%
%Choose EPOC
EPOC = input('Please enter the EPOC to filter on (start, Choice, Collect, Trial [shock test only]: \n','s');


%% Create table for WinStay/LoseShift data from all mice
%WinStayLoseShift=table(animalNames,
% WinStayLoseShiftTable=table(Wins,Losses,WinStay,LoseShift,WinStayPer,LoseShiftPer);
% Headers2 ={'Wins','Losses','WinStay','LoseShift','WinStayPer','LoseShiftPer'};
% WinStayLoseShiftTable.Properties.VariableNames([1:6])=Headers2;

%%
%run TDTbin2mat on all of these, fill in the table, and store the datastructures in workspace
%for each animal....
% EPOC = 'Choice';

zall_struct = struct;

zall_struct.data = {};

rows=0;
for aaa=1:numAnimals
    animalName=Table.animalNames{aaa};
    if contains(animalName, '-')
        animalName = strrep(animalName, '-', '_');
    end
    blockpath=Table.blockpath{aaa};
    whichStreams=Table.whichStreams(aaa);
    whichTTL=Table.whichTTL(aaa);
    behavFiles=Table.behavFiles{aaa};
    boris_files=Table.boris_files{aaa};
    SLEAP_files=Table.SLEAP_files{aaa};
    session_num = Table.session(aaa);
    SLEAP_time_range_adjustment=Table.SLEAP_time_range_adjustment{aaa};
    [data, meanZAnimal,zerrAnimal,times1,numTrials,Chan405,Chan465,TTLname, Descriptives, BehavData, zall, zall_motion, Y_dF_all, uv_temp, SLEAP_data_vel_filtered_session, Y_dF_all_session] =Analysis_BL2initiation_fn_v7(animalName,blockpath,whichStreams,whichTTL,behavFiles, EPOC, boris_files, SLEAP_files, SLEAP_time_range_adjustment, which_sessions);
    if which_sessions == '4'

    elseif which_sessions ~= '4'
            WinStayLoseShiftTable(aaa,:)=Descriptives;
    end
    if meanZAnimal==999
        
        disp(['Animal' Table.animalNames{aaa} 'has none of the specified trials'])
        Table.numofTrials(aaa)={numTrials};
    else
        ts1=times1;
          rows=rows+1;
        %change NEWNAME=zAnimalAll
        Table.Channel_405_name(aaa)={Chan405};
        Table.Channel_465_name(aaa)={Chan465};
        Table.TTL_name(aaa)={TTLname};
        Table.numofTrials(aaa)={numTrials};
%         WinStayLoseShiftTable(aaa,:)=Descriptives;

        %the correction below will adjust data points according to the first
        %inputed dataset.  This needs to be changed in case later sessions are
        %shorter (rather than longer) than the first.  THis error is only a
        %difference of one value - it's just because of the way the data are
        %truncated.
        if rows>1 
            if numel(meanZAnimal)> numel(zMultiAnimal(1,:))
                disp('time vector has one more element')
                newCol=numel(zMultiAnimal(1,:));
                meanZAnimal=meanZAnimal(1,1:newCol);
                ts1=times1(1,1:newCol);
                zerrAnimal=zerrAnimal(1,1:newCol);
            elseif numel(meanZAnimal)<numel(zMultiAnimal(1,:))
                zMultiAnimal=zMultiAnimal(:,1:numel(meanZAnimal));
                ZerrMultiAnimal= ZerrMultiAnimal(:,1:numel(meanZAnimal));
               
                disp('time vector has one less element')
            end
        end

        zMultiAnimal(rows,:)=meanZAnimal;
%         ZerrMultiAnimal(rows,:)=zerrAnimal;
%         master_struc(aaa) = data;
        behav_struct(aaa) = struct(BehavData);
        Y_dF_all_struct.data(aaa) = {Y_dF_all};

        session_num_field = append('Session', '_', string(session_num));
        uv.(animalName).(session_num_field) = uv_temp;
        zall_struct.data(aaa)={zall};
        zall_motion_struct.data(aaa)={zall_motion};
        ts1_struct.data = {ts1};
        Y_dF_all_session_struct.data(aaa) = {Y_dF_all_session};
        SLEAP_data_vel_session.data(aaa) = {SLEAP_data_vel_filtered_session};
    end
    
   
end
%% Concatenate IDs with WS/LS data

if which_sessions == '4'
    corrTable = Table;
elseif which_sessions ~= '4'
    WSLS_Table = [ID_Table, WinStayLoseShiftTable];
    NewTable = [Table, WinStayLoseShiftTable];
    %convert cell to matrix to be able to filter out animals with 0 trials
    NewTable.numofTrials = cell2mat(NewTable.numofTrials);
    searchval = 999; %numofTrials is set to 999 if an animal has 0 trials
    corrTable = NewTable(~any(NewTable.numofTrials == searchval, 2), :); %create a new table, corrTable, for correlations that only includes animals with trials
end

%%
%make user check table



%collect z scores

%%
[~,cols]=size(zMultiAnimal);

for ii=1:cols
    zAvg(ii)=mean(zMultiAnimal(:,ii));
    values=zMultiAnimal(:,ii);
    stdDev(ii)=std(values);
end

% for ii=1:cols
%     zAvgNAN(ii)=nanmean(zMultiAnimal(:,ii));
%     values=zMultiAnimal(:,ii);
%     stdDev(ii)=std(values);
% end
% 


SEM=stdDev./(sqrt(numAnimals));
lo=zAvg-SEM;
hi=zAvg+SEM;

%for TRANGE = -31 41
% ts1forprism = ts1(2137:3967);
% zAvgforprism = zAvg(2137:3967);
% semforprism = SEM(2137:3967);

%look for timeseries values > -10 but < 10, pull the relevant values from
%these variables to an array to plot in MATLAB or GraphPad Prism. Change
%the values from -10 / 10 to whatever is necessary
ind = ts1(1,:) > -10 & ts1(1,:) < 10;
ts1forprism = ts1(:,ind);
zAvgforprism = zAvg(:,ind);
semforprism = SEM(:,ind);
%for TRANGE = -41 51, gives you -10 to +10 around EPOC
% ts1forprism = ts1(3153:5188);
% zAvgforprism = zAvg(3153:5188);
% semforprism = SEM(3153:5188);

%FOR TRANGE = -51 61, gives you -8 to +8 around EPOC
% zAvgforprism_longer_window = zAvg(4374:6002);
% ts1forprism_longer_window = ts1(4374:6002);
% semforprism_longer_window = SEM(4374:6002);



%%
figure(1)
hold on



% Plot vertical line at epoch onset, time = 0

%for loop for 5 tone retrieval- retrieval is 50 tones though??

% a=0;
% for ii=1:50
%     p1=patch([a a+30 a+30 a], [-5 -5 20 20], [.8 1 1], 'EdgeColor','none');
%     a=a+35;
% end 

% vline = plot([0 0],[min(lo) max(hi)],'--','LineWidth',2)

% y = ylim; % current y-axis limits
vline = plot([0 0],[min(lo) max(hi)]);
vline.Color = 'k';
vline.LineStyle = '--';
vline.LineWidth = 2;

xxx=[ts1, ts1(end:-1:1)];
yyy=[lo, hi(end:-1:1)];
hp= fill(xxx,yyy,[ .8 .8 .8]); hold on;
set(hp,'FaceColor', [ .8 .8 .8],'EdgeColor','none','FaceAlpha',0.8,'EdgeAlpha',0.8);hold on;
p3=plot(ts1, zAvg, 'k', 'LineWidth', 2,'DisplayName','z-score'); hold on;



% Make a legend
% legend([p3 hp],{'Block 1','SEM'}, 'AutoUpdate', 'off');

set(gca,'Layer','top')

%%
c = distinguishable_colors(numel(ID_Table));

figure(2);
set(gca, 'ColorOrder', c, 'NextPlot', 'replacechildren');
% % Create 10 data series
% XY = rand(1000,10);
% % Build colormap and shuffle
% cmap = colormap(jet(size(XY,2)));
% cmap = cmap(randperm(length(cmap)),:)
% %Set colororder and plot
% ax = axes('colororder',cmap);hold on






plot(ts1, zMultiAnimal, 'LineWidth',2);

ID_legend = table2cell(ID_Table);
legend(ID_legend);

%% Truncate zMultiAnimal to just the window we care about quantifying: better way to do this below
% ind2 = ts1(1,:) > -4 & ts1(1,:) < 0;
% zMultiAnimal_Ind = zMultiAnimal(:,ind2);
% for ii = 1:size(zMultiAnimal_Ind,1);
%     zTrunc(ii,1) = mean(zMultiAnimal_Ind(ii)); 
%     ii = ii+1;
% end

%% Calculate AUC for different time windows, store in AUC array
AUC=[]; % cue, shock
for qq = 1:size(zMultiAnimal,1);
    AUC(qq,1)=trapz(zMultiAnimal(qq,ts1(1,:) < -0 & ts1(1,:) > -2));
    AUC(qq,2)=trapz(zMultiAnimal(qq,ts1(1,:) > 0 & ts1(1,:) < 5));
    qq=qq+1;
end
AUC_mean = mean(AUC);
AUC_std = std(AUC)/sqrt(numel(AUC(:,1)));
%%
means=[]; % cue, shock
for qq = 1:size(zMultiAnimal,1);
    means(qq,1)=mean(zMultiAnimal(qq,ts1(1,:) < -0 & ts1(1,:) > -2));
    means(qq,2)=mean(zMultiAnimal(qq,ts1(1,:) > 0 & ts1(1,:) < 5));
    qq=qq+1;
end
means_mean = mean(means);
means_std = std(means)/sqrt(numel(means(:,1)));


%%
%if a given animal has 0 trials, will need to delete rows to deal with
%mismatch
figure;
scatter(corrTable.RiskPercent,AUC(:,1));
coefficients = polyfit(corrTable.RiskPercent, AUC(:,1), 1);
xFit = linspace(0, max(corrTable.RiskPercent), 1000); %min(WSLS_Table.RiskPercent)
yFit = polyval(coefficients , xFit);
hold on;
title('Scatter: Risk % vs Pre-Choice')
xlabel('Proportion Risky Choices')
ylabel('Normalized Pre-Choice AUC')
plot(xFit, yFit, 'r-', 'LineWidth', 2);
grid on;

figure;
scatter(corrTable.RiskPercent,AUC(:,2));
coefficients = polyfit(corrTable.RiskPercent, AUC(:,2), 1);
xFit = linspace(0, max(corrTable.RiskPercent), 1000); %min(WSLS_Table.RiskPercent)
yFit = polyval(coefficients , xFit);
hold on;
title('Scatter: Risk % vs Post-Choice')
xlabel('Proportion Risky Choices')
ylabel('Normalized Post-Choice AUC')
plot(xFit, yFit, 'r-', 'LineWidth', 2);
grid on;

figure;
scatter(corrTable.LoseShiftPercent,AUC(:,1));
coefficients = polyfit(corrTable.LoseShiftPercent, AUC(:,1), 1);
xFit = linspace(0, max(corrTable.LoseShiftPercent), 1000); %min(WSLS_Table.RiskPercent)
yFit = polyval(coefficients , xFit);
hold on;
title('Scatter: LoseShift vs Pre-Choice')
xlabel('Proportion of Lose Shifts')
ylabel('Normalized Pre-Choice AUC')
plot(xFit, yFit, 'r-', 'LineWidth', 2);
grid on;

%%
figure;
C = [.8 .2 .8;
     .2 .75 .2];
CE = [.5 .1 .5];
superbar(AUC_mean, 'BarFaceColor', C, 'BarWidth', .5, 'E', AUC_std, 'ErrorbarColor', CE);
xlim([0.5 2.5]);

% title('superbar');

%% IF USING ALL SESSIONS + Need to double check that any mice with 1 session only are removed

LoseShift_D(:,1) = corrTable.LoseShiftPercent(1:2:end);
LoseShift_D(:,2) = corrTable.LoseShiftPercent(2:2:end);
WinStay_D(:,1) = corrTable.WinStayPercent(1:2:end);
WinStay_D(:,2) = corrTable.WinStayPercent(2:2:end);
RiskyPercent_D(:,1) = corrTable.RiskPercent(1:2:end);
RiskyPercent_D(:,2) = corrTable.RiskPercent(2:2:end);
meanRiskyPercent_D = mean(RiskyPercent_D);
SEMRiskyPercent_D = std(RiskyPercent_D)/sqrt(numel(RiskyPercent_D(:,1)));

%%
figure;
C = [.8 .2 .8;
     .2 .75 .2];
CE = [.5 .1 .5];
superbar(meanRiskyPercent_D, 'BarFaceColor', C, 'BarWidth', .5, 'E', SEMRiskyPercent_D, 'ErrorbarColor', CE);
xlim([0.5 2.5]);
% title('superbar');
%% testing from Reed 2018
% findchangepts(zMultiAnimal(:,ts1(1,:) > 0 & ts1(1,:) < 2))

%%
% Specify the variables you want to save
variables_to_save = {
    'SLEAP_data_vel_session', 
    'Y_dF_all_session_struct', 
    'Y_dF_all_struct', 
    'behav_struct', 
    'corrTable', 
    'ts1_struct', 
    'uv', 
    'zall_motion_struct', 
    'zall_struct'
};

% Open a dialog box to select the save location
[file, path] = uiputfile('*.mat', 'Save Workspace Variables');

% Check if the user canceled the dialog
if isequal(file,0) || isequal(path,0)
    disp('User canceled the operation');
else
    % Construct the full file path
    full_path = fullfile(path, file);
    
    % Create a structure containing the variables you want to save
    save_data = struct();
    for i = 1:numel(variables_to_save)
        var_name = variables_to_save{i};
        if exist(var_name, 'var')
            save_data.(var_name) = eval(var_name);
        end
    end
    
    % Save the selected variables to the specified .MAT file
    save(full_path, '-struct', 'save_data');
    disp(['Saved workspace variables to ' full_path]);
end