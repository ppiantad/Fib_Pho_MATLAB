[BehavData,ABETfile,Descriptives]=ABET2TableFn_Chamber_Av3('159 12182020.csv');
[~,~,boris] = xlsread('68 RDT 20191031_checked.csv');
boris_extract = boris(17:end,[1,6]);
boris_Extract_tbl = cell2table(boris_extract)
boris_Extract_tbl.Properties.VariableNames = {'choiceTime','type'};
outerjoin(boris_Extract_tbl, BehavData);

abet_and_boris = tblvertcat(boris_Extract_tbl, BehavData);
tblB = sortrows(abet_and_boris,{'choiceTime'},{'ascend'});