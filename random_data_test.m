C19 = readmatrix("E:\MATLAB\TDTbin2mat\Inscopix Python\BLA-Insc-3\SingleCellAlignmentData\C16\Block_Reward Size_Choice Time (s)\(2.0, 'Large')\plot_ready.csv");
zscore = [-10 -5]
window = [-10 10];
ts1 = window(1):0.1:window(2)-0.1;
C19_2 = C19(2:end,2:end);

zall = zeros(size(C19_2));
tmp = 0;
pp = 1;


for i = 1:size(C19_2,1)
%     BL_shifted(pp,:)=[BASELINE_PER(1)+time2EPOC(i) BASELINE_PER(2)+time2EPOC(i)]; %BL_shifted(pp,:)=[BASELINE_PER(1)+(-1*time2Collect(i)) BASELINE_PER(2)+time2Collect(i)];
%     ind = ts2(1,:) < BL_shifted(pp,2) & ts2(1,:) > BL_shifted(pp,1);
    bl_ind = ts1(1,:) <= zscore(2) & ts1(1,:) >= zscore(1);
    
    %use if you want to take the Z-score using the entire window mean
%     zb = mean(C19_2(i,:)); % baseline period mean
%     zsd = std(C19_2(i,:)); % baseline period stdev
%     
    %use if you want to calculate the Z-score using your specified baseline
    zb = mean(C19_2(i,bl_ind)); % baseline period mean
    zbmedian = median(C19_2(i,length(ts1)));
    zsd = std(C19_2(i,bl_ind)); % baseline period stdev




    pp=pp+1;
    for j = 1:size(C19_2,2) % Z score per bin
        tmp = tmp + 1;
        zall(i,tmp)=(C19_2(i,j) - zb)/zsd;
%         dfALL(i,tmp)=(Y_dF_all(i,j) - zbmedian)/zbmedian;
        BL_mean(i) = mean(C19_2(i,:));
    end
    tmp=0;
end


zerror = std(zall)/size(zall,1);
ZallMean=mean(zall,1);