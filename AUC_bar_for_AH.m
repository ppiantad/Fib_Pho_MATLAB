
f = figure;  
f.Position = [10 10 200 300]; 

bar([AUC_mean(2,2), AUC_mean(2,1)]);
hold on;
er = errorbar([AUC_mean(2,2), AUC_mean(2,1)],[AUC_sem(2,2),AUC_sem(2,1)]);    
er.Color = [0 0 0];    
fontsize(f, 14, "points")
set(gca,'XTick',[])
er.LineStyle = 'none';  
hold off;

% get stats from the AUCs that you created. as of 5/8/2023 (:,1) should compare Pre-Choice
% AUC
[h,p,ci,stats] = ttest2(AUC_iters{1, 1}(:,1),AUC_iters{1, 2}(:,1));


%%

f = figure;  
f.Position = [10 10 200 300]; 

bar([AUC_mean(1,2), AUC_mean(1,1)]);
hold on;
er = errorbar([AUC_mean(1,2), AUC_mean(1,1)],[AUC_sem(1,2),AUC_sem(1,1)]);    
er.Color = [0 0 0];  
fontsize(f, 14, "points")
set(gca,'XTick',[])
er.LineStyle = 'none';  
hold off;

% get stats from the AUCs that you created. as of 5/8/2023 (:,2) should compare Post-Choice
% AUC
[h,p,ci,stats] = ttest2(AUC_iters{1, 1}(:,2),AUC_iters{1, 2}(:,2));

format long
p_decimal = fprintf('%.4f\n', p);