cd('E:\MATLAB\TDTbin2mat\Inscopix Python\Individual_Mouse_Cell_Data')

folder = cd;
FileList = dir(fullfile(folder,'**'));

inscopix_risk_struct = struct;
name_substring = "BLA-Insc";

jj = 1;
match = '-';
new = '_';

for ii = 1:size(FileList, 1)
 if contains(FileList(ii).name,name_substring,'IgnoreCase',true) == true
     IDs{jj,:} = FileList(ii).name;
     Mouse_IDs{jj,:} = replace(IDs{jj}, match, new);
     inscopix_risk_struct.(Mouse_IDs{jj}) = struct;
     jj = jj+1;
 end
    
end

%go to folder w/ Inscopix Data
topLevelFolder = pwd;
files = dir(topLevelFolder);
dirFlags = [files.isdir];
subFolders = files(dirFlags);
subFolderNames = {subFolders(3:end).name}
subFolderNames_regexp = regexprep(subFolderNames, '\W','_');

for ii = 1:size(Mouse_IDs, 1)
    inscopix_risk_struct.(Mouse_IDs{ii}).SingleCellAlignmentData = struct;

end


%navigate to directory w/ list of Cells @ run topLevelFolder above in each
%mouse's SingleCellAlignmentData folder
for ii = 1:size(subFolderNames, 2)
    inscopix_risk_struct.BLA_Insc_1.SingleCellAlignmentData.(subFolderNames{1,ii}) = struct;
end

for ii = 1:size(Mouse_IDs, 1)
    for jj = 1:size(fieldnames(inscopix_risk_struct.(Mouse_IDs{ii}).SingleCellAlignmentData),1)
        fieldnames_list = fieldnames(inscopix_risk_struct.(Mouse_IDs{ii}).SingleCellAlignmentData);
        for kk = 1:size(subFolderNames, 2)
            inscopix_risk_struct.(Mouse_IDs{ii}).SingleCellAlignmentData.(fieldnames_list{jj}).(regexprep(subFolderNames{1,kk}, '\W','_')) = struct;
            
        end
    end
end

for ii = 1:size(fieldnames(inscopix_risk_struct.BLA_Insc_1.SingleCellAlignmentData.C01),1)
    
end

%% https://www.mathworks.com/matlabcentral/answers/522194-how-to-create-a-loop-that-runs-a-function-and-save-output-through-subfolders-in-a-directory
parentDir = pwd;
FolderStructure = dir;
con = struct2cell(FolderStructure);
myFiles = con(1,:)';
Phases_I=[];
Phases_O=[];
GR_I=[];
GR_O=[];

for ii = 1:size(Mouse_IDs, 1)
    for k = 1:size(fieldnames(inscopix_risk_struct.BLA_Insc_1.SingleCellAlignmentData.C01),1)
        for z = 3:length(myFiles)
            name = string(myFiles(z));
            dr = ['E:\MATLAB\TDTbin2mat\Inscopix Python\Individual_Mouse_Cell_Data\BLA-Insc-1\SingleCellAlignmentData\C01\',char(name)];
            fs=dir(dr);
            fs=fs(~ismember({fs.name},{'.','..'})); %dir adds 2 files '.' and '..' that correspond to the current and previous folder. remove these
            cd(dr);
            fs=struct2cell(fs);
            subF= fs(1,:)';
            
            
                for p = 1:length(subF)
                    if contains(subF{p},'.0') == true
                        subF{p} = regexprep(subF{p}, '\,','');
                        subF{p}=regexprep(subF{p}, '\,','');
                        subF{p}=regexprep(subF{p}, '\,','');
                        subF{p}=regexprep(subF{p}, '\(','');
                        subF{p}=regexprep(subF{p}, '\)','');
                        subF{p}=regexprep(subF{p}, '\''','');
                        subF{p}=regexprep(subF{p}, '\.0','');
                        subF{p}= regexprep(subF{p}, '\W','_');
                        subF{p} = append('Block_',subF{p});
                        inscopix_risk_struct.BLA_Insc_1.SingleCellAlignmentData.C01.(subFolderNames_regexp{z}).(subF{p}) = struct;
                    elseif contains(subF{p},'.0') == false
                        inscopix_risk_struct.BLA_Insc_1.SingleCellAlignmentData.C01.(subFolderNames_regexp{z}).(subF{p}) = struct;
                    end 
                        
                    
                end
        end
         
                
    end
            
%             for k = 3:length(subF1)
%                 
%                 
%             end
%             
            
        end
end
%% this allows you to add new levels etc to your data structure (for example, adding a meta structure called mouseID. modify as necessary
inscopix_risk_struct_added = struct('mouseID',struct(inscopix_risk_struct));

