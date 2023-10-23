%% set root directory (where files are located on PC)
root = upper('d:');

%%

% behavFiles = dir('D:\MATLAB\TDTbin2mat\Radke_test'); 
% 
% for i = 1:length(myFiles) % for each file
%     file = myFiles(i) % get the file name
%     behavFiles(i)=myFiles.name{i}
%     i=i+1;
%     % insert your code
%     
% end
% 
% 
% behavFiles={
%     'LAPTOP-GQB8OREV_HTR2B_Reversal (Marble) vHolmes (Mouse Pairwise Discrimination v3)_304.csv'...
%     'LAPTOP-GQB8OREV_HTR2B_Reversal (Marble) vHolmes (Mouse Pairwise Discrimination v3)_321.csv'...
%     };
% 
% numFiles = numel(behavFiles);

%%
xlfiles = dir('*.csv'); % Make sure "Current Folder" is where you have all the .csvs that you want to analyze 
Nfiles = length(xlfiles) ;  % number of files
% loop for each file 
for i = 1:Nfiles
    fname = xlfiles(i).name ;   % file name
    data = xlsread(fname) ;     % read the file 
    behavFiles = {xlfiles.name}; %.'
    %%do what you want %%%
end
Genotypes = readtable('Genotypes.txt');
numFiles = numel(behavFiles);
%%
MouseID=behavFiles;
GroupID=behavFiles;
DateTime=behavFiles;
SessionID=behavFiles;
TotalWins=zeros(1,numFiles);
TotalLosses = zeros(1,numFiles);
TotalWinStay = zeros(1,numFiles);
TotalLoseShift =zeros(1,numFiles);
WinStayPercent = zeros(1,numFiles);
LoseShiftPercent = zeros(1,numFiles);
% ID_Table = cell2table(behavFiles);

%%
%make table with all animal info
Table=table(behavFiles',MouseID',GroupID',DateTime',SessionID',TotalWins',TotalLosses',TotalWinStay',TotalLoseShift',WinStayPercent',LoseShiftPercent');
Headers={'behavFiles','MouseID','GroupID','DateTime','SessionID','TotalWins','TotalLosses','TotalWinStay','TotalLoseShift','WinStayPercent','LoseShiftPercent'};
Table.Properties.VariableNames([1:11])=Headers;


%%
rows=0;
for aaa=1:numFiles
    behavFiles=Table.behavFiles{aaa};
    [data, ABETdata, Descriptives] = ABET2TableFn_Radke_fn_v1(behavFiles, Genotypes);
        
        Table.MouseID(aaa)={Descriptives.MouseID};
        Table.GroupID(aaa)={Descriptives.GroupID};
        Table.DateTime(aaa)={Descriptives.DateTime};
        Table.SessionID(aaa)={Descriptives.SessionID};
%         Table(aaa,(5:10)) = (Descriptives(aaa,5:10));
        Table.TotalWins(aaa)=(Descriptives.TotalWins);
        Table.TotalLosses(aaa)=(Descriptives.TotalLosses);
        Table.TotalWinStay(aaa)=(Descriptives.TotalWinStay);
        Table.TotalLoseShift(aaa)=(Descriptives.TotalLoseShift);
        Table.WinStayPercent(aaa)=(Descriptives.WinStayPercent);
        Table.LoseShiftPercent(aaa)=(Descriptives.LoseShiftPercent);
        Table.Genotype(aaa)=(Descriptives.Genotype);
%         WinStayLoseShiftTable(aaa, 1)=Descriptives;

    rows=rows+1;
end


%% Sort by MouseID
Table = sortrows(Table,'MouseID','ascend');
%%
%random stuff to try to organize table by date
% datetimesplit = regexp(Table.DateTime,' ', 'split');
% datetimeextracted = arrayfun(@(n) datetimesplit{n}(1),1:size(Table.DateTime,1))' %//'
% datetimeextracted = datetime(datetimeextracted,'InputFormat','MM/dd/yyyy');
% Table.DateTime2 = datetimeextracted;
% Table = sortrows(Table,'DateTime2','ascend');
%%
%Early, Mid, Late
% early = floor(rows/3);
% mid = ceil(rows/3);
% late = floor(rows/3);
% 
% for aaa=1:early
%     if WinStayLoseShiftTable.MouseID(aaa)
%         
%     end
% end


%%
% Uncomment "writetable" to save WSLS Table as a .xlsx
% writetable(Table,'WSLS_all.xlsx')
   