function [data, trials, varargin, mask, k] = TrialFilter_v2(data,varargin);
VALID_PARS = {'BLOCK','TYPE','SHK','REW','OMIT','OMITALL','ALL','WSLS','STTOCHO','WINSTAY','LOSESHIFT','LOSEOMIT', 'WIN', 'LOSS'};
%for varargin
%JSFilter will assign the valid parameter equal to the argument that
%immediately succeeds that paramter.  So, for example, if you call the
%function as follows: d=JSFilter(data,BigRew_choice,'SESSION',SESSION3), it
%will set SESSION=Session3.  Parameters should be entered as strings, and
%assignment values followng the parameters should be vectors.


%initialize variables
BLOCK=[]; %BLOCK number (1,2,3,4,5)
TYPE=[]; %use 1 for force, 0 for free
SHK=[]; %use 0 for no shock, 1 for shock
REW=[]; % use 1.2 for big, 0.3 for small
OMIT=[]; %use 0 for non-omissions, 1 for free trial omissions, 2 for force trial omissions
ALL=[]; %if you want to see all trials, set ALL to 1 (or any number should work)
WSLS=[]; %1=win; 2=win+1 (trial after win); 3=loss; 4=loss+1;
WINSTAY=[];
LOSESHIFT=[];
LOSEOMIT=[];
WIN=[]; % 1 = win
LOSS=[]; % 3 = loss
AA=[];



% parse varargin
for ii = 1:2:length(varargin)
    if ~ismember(upper(varargin{ii}), VALID_PARS)
        error('%s is not a valid parameter', upper(varargin{ii}));
    end
    eval([upper(varargin{ii}) '=varargin{ii+1};']);  %sets the argument in equal to the value that follows it.  

end




names = {'one','two','three','four'}
for index = 1: length(names)
    mask2.(names{index}) = index;
end



for k=1:2:length(varargin)
    %Filtering by SESSION: SESSION will be set equal to the vector which
    %corresponds to the SESSION of interest.
    if strcmp(upper(varargin{k}),'BLOCK')
        block_label = [varargin{k} num2str(varargin{k+1})];
        trials{k} = table2cell(data([find(data.Block==BLOCK)],1));
        mask.(block_label) = ismember(data.Trial,cell2mat(trials{1,k}))
%         data.(varargin{k}) = data.Block==BLOCK;
%         trials_test = ismember(BehavData.Trial,cell2mat(trials{1,1}))
%         trials3 = {table2cell(BehavData([find(BehavData.Block==1)],1))};

    end
    
    if strcmp(upper(varargin{k}),'TYPE')
        trials{k} = table2cell(data([find(data.ForceFree==TYPE)],1));
        mask.(varargin{k}) = ismember(data.Trial,cell2mat(trials{1,k}))
%         data([find(data.ForceFree~=TYPE)],:)=[];
    end
    
    if strcmp(upper(varargin{k}),'SHK')
        trials{k} = table2cell(data([find(data.shock==SHK)],1));
        mask.(varargin{k}) = ismember(data.Trial,cell2mat(trials{1,k}))
%         data([find(data.shock~=SHK)],:)=[];
    end
    
    
    if strcmp(upper(varargin{k}),'REW')
        
        trials{k} = table2cell(data([find(data.bigSmall==REW)],1));
        mask.(varargin{k}) = ismember(data.Trial,cell2mat(trials{1,k}))
        
%        trials = table2cell(data([find(data.bigSmall==REW)],1));
%        data.(varargin{k}) = data.bigSmall==REW;
%        data([find(data.bigSmall~=REW)],:)=[];
    end
    
    if strcmp(upper(varargin{k}),'OMIT')
        trials{k} = table2cell(data([find(data.omission==OMIT)],1));
        mask.(varargin{k}) = ismember(data.Trial,cell2mat(trials{1,k}))
%         data([find(data.omission~=OMIT)],:)=[];
    end
    
    if strcmp(upper(varargin{k}),'OMITALL')
        trials{k} = table2cell(data([find(data.omissionALL==OMITALL)],1));
        mask.(varargin{k}) = ismember(data.Trial,cell2mat(trials{1,k}))
%         data([find(data.omissionALL~=OMITALL)],:)=[];
    end
    
    if strcmp(upper(varargin{k}),'WSLS')
        trials{k} = table2cell(data([find(data.WSLScode==WSLS)],1));
        mask.(varargin{k}) = ismember(data.Trial,cell2mat(trials{1,k}))
%         data([find(data.WSLScode~=WSLS)],:)=[];
    end
    
    if strcmp(upper(varargin{k}),'ALL')
       trials{k} = num2cell(data.Trial);
       mask.(varargin{k}) = ismember(data.Trial,cell2mat(trials{1,k}))
       disp('Using all trials')
    end
    
    if strcmp(upper(varargin{k}),'WINSTAY')
        trials{k} = table2cell(data([find(data.win_stay==WINSTAY)],1));
        mask.(varargin{k}) = ismember(data.Trial,cell2mat(trials{1,k}))
%         data([find(data.win_stay~=WINSTAY)],:)=[];
    end
    
    if strcmp(upper(varargin{k}),'LOSESHIFT')
        trials = table2cell(data([find(data.lose_shift==LOSESHIFT)],1));
        mask.(varargin{k}) = ismember(data.Trial,cell2mat(trials{1,k}))
%         data([find(data.lose_shift~=LOSESHIFT)],:)=[];
    end
    
    if strcmp(upper(varargin{k}),'LOSEOMIT')
        trials{k} = table2cell(data([find(data.lose_omit==LOSEOMIT)],1));
        mask.(varargin{k}) = ismember(data.Trial,cell2mat(trials{1,k}))
%         data([find(data.lose_omit~=LOSEOMIT)],:)=[];
    end    
    
    if strcmp(upper(varargin{k}),'WIN')
        trials{k} = table2cell(data([find(data.WL==WIN)],1));
        mask.(varargin{k}) = ismember(data.Trial,cell2mat(trials{1,k}))
%         data([find(data.WL~=WIN)],:)=[];
    end    
    
    if strcmp(upper(varargin{k}),'LOSS')
        trials{k} = table2cell(data([find(data.WL==LOSS)],1));
        mask.(varargin{k}) = ismember(data.Trial,cell2mat(trials{1,k}))
%         data([find(data.WL~=LOSS)],:)=[];
    end   
    
    if strcmp(upper(varargin{k}),'LOSS')
        trials{k} = table2cell(data([find(data.WL==LOSS)],1));
        mask.(varargin{k}) = ismember(data.Trial,cell2mat(trials{1,k}))
%         data([find(data.WL~=LOSS)],:)=[];
    end   
end

iter = 1;
filter_tbl = struct2table(mask);

for ii=1:size(filter_tbl,1);
    if filter_tbl{ii,iter} == 1
        if filter_tbl.
            whichBlock=regexp(ABETdata{ii,4},'n','split'); %pull block from "SX-Free...=big" or similar name
            currBlock=str2num(whichBlock{2});


out1 = filter_tbl.BLOCK1 == 1 & filter_tbl.REW == 1
out2 = out1 == 1 & filter_tbl.TYPE == 1
if filter_tbl.BLOCK1(1) == 1;
    & filter.

% for yy = 1:size(filter_tbl,1);
%     
%     trials = ismember(filter_tbl(yy,:)
    
% results_full = nan(size(data,1));
% 
% [temp,indx_full,indx_order] = intersect(mask2.BLOCK,mask2.REW);
% results_full(indx_full) = mask2.REW(indx_order);
% results_full
% trials = ismember(cell2mat(trials{1,1}),cell2mat(trials{1,3}))
% data([find(data.Block~=BLOCK)],:)=[];
% data([find(data.ForceFree~=TYPE)],:)=[];
%         mask2 = data.Trial == cell2mat(trials);
%         data.FilterCol(mask2) = 1;
end

