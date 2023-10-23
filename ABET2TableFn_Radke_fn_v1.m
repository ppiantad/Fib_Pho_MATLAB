function [data, ABETdata, Descriptives] = ABET2TableFn_Radke_fn_v1(behavFiles, Genotypes)


%ABET2Table creates a table with columns (outlined under "column headers"
%comment) from ABET behavioral data "behavFiles"

%ABET file should be either "Raw Data" extracted from ABET. This is
%currently written such that data comes from the "batch" version of the
%ABET data extraction, which adds 18 lines of headers to the file (with
%identifying info about the subject etc.

%column headers
% (1)Trial: trial number
% (2)Result: 1,2,3,4, or 5
% (3)stTime: trial start time
% (4)choicetime: timestamp of choice 
% (5)collection time: timestamp of reward collection
% (10)WL: 1 if win; 3 if loss
% (11)WSLScode: 2 if win+1; 4 if loss+1;

%initiate table
[~,~,ABETdata]=xlsread(behavFiles);

Headers={'Trial','Result','stTime','choiceTime','collectionTime','WL','WSLScode','win_stay','lose_shift'};
data=table(zeros(80,1),zeros(80,1),zeros(80,1),zeros(80,1),zeros(80,1),zeros(80,1),zeros(80,1), zeros(80,1), zeros(80,1));
data.Properties.VariableNames([1:9])=Headers;

%%
%add ABET data to table
%loop through all rows
[rows,~]=size(ABETdata);
trial=1;

%find the first nonzero timestamp (all timestamps at 0 are program checks,
%and we don't care about these when we're searching for behavioral events
%rr = 18 because the files pulled in batch by the most recent version of
%ABET contain a bunch of headers that are irrelevant, don't need to scan
%these. this loop is probably unnecessary (could just set rr = 29), but
%will keep it this way because could be useful for future modification. 

stop=999; rr=18;
while stop>0
    startRow=rr;
    if ABETdata{rr,2} == 29

        stop=-999;
    end
    rr=rr+1;
    
end

%%
% string names for all possible trial outcomes
corr_str = 'Correct';
inc_str = 'Incorrect';
inc_corr_str = 'Incorrect Correction Trial';
corr_corr_str = 'Correct Correction Trial';

%loop through all rows of the ABET file, extracting the relevant timestamps
%and labeling each trial Result

for ii=startRow:rows
    
    data.Trial(trial)=trial;
    
    %keep track of Result
    % 1 = correct trial
    % 2 = incorrect trial
    % 3 = incorrect correction trial
    % 4 = correct correction trial
    if strcmp(ABETdata{ii,4}, corr_str)
        data.Result(trial) = 1;
    end
    
    if strcmp(ABETdata{ii,4}, inc_str)
        data.Result(trial) = 2;
    end
    
    if strcmp(ABETdata{ii,4}, inc_corr_str)
        data.Result(trial) = 3;
    end

    if strcmp(ABETdata{ii,4}, corr_corr_str)
        data.Result(trial) = 4;
    end
    
    %TrialStart
    if strcmp(ABETdata{ii,4},'Next trial')
        data.stTime(trial)=ABETdata{ii,1};
    end
    
     %CHOICE TIME
    %if it's a choice
    if ABETdata{ii,2}==1 && ABETdata{ii,6} == 7
        data.choiceTime(trial)=ABETdata{ii,1};  
        trial = trial+1;
    end
   
    
    %COLLECTION TIME
    %because I increment the trial based on choiceTime, add each
    %reward retrieved to the previous trial
    if regexp(ABETdata{ii,4},'Reward Collected*')
        data.collectionTime(trial - 1)=ABETdata{ii,1};
    end
    
end

%add win stay/lose shift info.  To do this, add a column for
%win-stay/lose-shift code, called WSLS code.  For this code,
% if trial is a win, code=1;
% if trial is the trial after a win, code=2;
% if trial is a loss, code=3;
% if trial is a trial after a loss, code=4; currenly this code counts
% correct correction trials as wins, and incorrect correction trials as
% losses, does not differentiate. possible to change later


for jj=1: numel(data.Trial)

    if data.Result(jj) == 1 || data.Result(jj) == 4
                data.WL(jj)=1; %win
            elseif data.Result(jj)== 2 || data.Result(jj)== 3
                data.WL(jj)=3; %loss
            end
%       if data.Result(jj) == 1
%                 data.WL(jj)=1; %win
%             elseif data.Result(jj)== 2 
%                 data.WL(jj)=3; %loss
%             elseif data.Result(jj) > 2
%                 data.WL(jj)=0;
%             end    
            
        if jj>1
            if data.WL(jj-1)==1
                data.WSLScode(jj)=2; %win+1 trial
                if data.Result(jj) == 1
                    data.win_stay(jj) = 1; %win_stay is 1 if chose big after a win
                end
                
            elseif data.WL(jj-1) == 2 || data.WL(jj-1) == 3
                data.WSLScode(jj) = 4; %loss+1 trial
               if data.Result(jj)== 4
                   data.lose_shift(jj)=1;
               end
               
            end
        end
  
end    


% deletes excess rows (if the default 80 rows are not filled); make sure
% this doesn't cause a problem on reversal days where there are a large #
% of trials
toDelete = data.Trial == 0;
data(toDelete,:) = [];
size(data);

%probably not necessary to create the tables this way, but it works
MouseID = upper(ABETdata{11,2});
GroupID = ABETdata{12,2};
DateTime = ABETdata{4,2};
SessionID = ABETdata{15,2};
TotalWins = sum(data.WL(:)==1);
TotalLosses = sum(data.WL(:)==3);
TotalWinStay = sum(data.win_stay(:)==1);
TotalLoseShift =sum(data.lose_shift(:)==1);
WinStayPercent = TotalWinStay / TotalWins;
LoseShiftPercent = TotalLoseShift / TotalLosses;

Descriptives = table;
Descriptives.MouseID = MouseID;
Descriptives.GroupID = GroupID;
Descriptives.DateTime = DateTime;
Descriptives.SessionID = SessionID;
Descriptives.TotalWins = TotalWins;
Descriptives.TotalLosses = TotalLosses;
Descriptives.TotalWinStay = TotalWinStay;
Descriptives.TotalLoseShift = TotalLoseShift;
Descriptives.WinStayPercent = WinStayPercent;
Descriptives.LoseShiftPercent = LoseShiftPercent;

% create a cell array from the Genotype table
Genotypes1 = table2cell(Genotypes(:,1:2));

% Loop through the table and compare the MouseID from the current loop to
% the Genotype cell array, pull the genotype from the 2nd column when the
% IDs match

for pp = 1:numel(Genotypes(:,1))
    if strcmp(Descriptives.MouseID, Genotypes1(pp,1));
       Descriptives.Genotype = Genotypes1(pp,2);
    end
end    



%%
% uncomment first one if you want to write descriptive statistics to file,
% unnecessary for now
% writetable(Descriptives,'filename.xlsx')







