%% set root directory (where files are located on PC)
root = upper('h:');

which_region = input('Which dataset do you want to analyze?\n 1) Code Test \n 2) BLA \n 3) vmPFC \n 4) D1 \n 5) D2 \n 6) BLA GFP \n 7) vmPFC GFP \n 8) GFP ALL \n 9) D1-iOP \n 10) D2-iOP \n 15) D1-eOP \n 16) D2-eOP \n 12) vHPC \n','s');

which_sessions = input('Which sessions do you want to analyze from these mice? \n 1) Early RM \n 2) Late RM \n 3) RDT \n','s');


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
rows=0;
for aaa=1:numAnimals
    animalName=Table.animalNames{aaa};
    blockpath=Table.blockpath{aaa};
    whichStreams=Table.whichStreams(aaa);
    whichTTL=Table.whichTTL(aaa);
    behavFiles=Table.behavFiles{aaa};
    boris_files=Table.boris_files{aaa};
    SLEAP_files=Table.SLEAP_files{aaa};
    SLEAP_time_range_adjustment=Table.SLEAP_time_range_adjustment{aaa};
    [timedownsample, Y_dF_all, block_end] = func_Photometry_Signal_Check_v3(blockpath, behavFiles, whichStreams, whichTTL);
    [block1_pks_per_min_filtered, block2_pks_per_min_filtered, block3_pks_per_min_filtered] = findpeaks_photometry_v1(timedownsample, Y_dF_all, block_end, animalName);
    pks_per_min_filtered(aaa,:) = [block1_pks_per_min_filtered block2_pks_per_min_filtered block3_pks_per_min_filtered];
%     WinStayLoseShiftTable(aaa,:)=Descriptives;
%     if meanZAnimal==999
%         
%         disp(['Animal' Table.animalNames{aaa} 'has none of the specified trials'])
%         Table.numofTrials(aaa)={numTrials};
%     else
%         ts1=times1;
%           rows=rows+1;
%         %change NEWNAME=zAnimalAll
%         Table.Channel_405_name(aaa)={Chan405};
%         Table.Channel_465_name(aaa)={Chan465};
%         Table.TTL_name(aaa)={TTLname};
%         Table.numofTrials(aaa)={numTrials};
% %         WinStayLoseShiftTable(aaa,:)=Descriptives;
% 
%         %the correction below will adjust data points according to the first
%         %inputed dataset.  This needs to be changed in case later sessions are
%         %shorter (rather than longer) than the first.  THis error is only a
%         %difference of one value - it's just because of the way the data are
%         %truncated.
%         if rows>1 
%             if numel(meanZAnimal)> numel(zMultiAnimal(1,:))
%                 disp('time vector has one more element')
%                 newCol=numel(zMultiAnimal(1,:));
%                 meanZAnimal=meanZAnimal(1,1:newCol);
%                 ts1=times1(1,1:newCol);
%                 zerrAnimal=zerrAnimal(1,1:newCol);
%             elseif numel(meanZAnimal)<numel(zMultiAnimal(1,:))
%                 zMultiAnimal=zMultiAnimal(:,1:numel(meanZAnimal));
%                 ZerrMultiAnimal= ZerrMultiAnimal(:,1:numel(meanZAnimal));
%                
%                 disp('time vector has one less element')
%             end
%         end
% 
%         zMultiAnimal(rows,:)=meanZAnimal;
% %         ZerrMultiAnimal(rows,:)=zerrAnimal;
% %         master_struc(aaa) = data;
%         behav_struct(aaa) = struct(BehavData);
%         Y_dF_all_struct.data(aaa) = {Y_dF_all};
% 
%         zall_struct.data(aaa)={zall};
%         zall_motion_struct.data(aaa)={zall_motion};
%         ts1_struct.data = {ts1};
%     end
%     
end