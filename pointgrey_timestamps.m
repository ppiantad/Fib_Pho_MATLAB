[NUM,TXT,RAW] = xlsread('102_RDT_D12020-01-15T13_09_09.csv');

TXT2 = TXT(:,13);

%%
b = cell(size(TXT2));
for i= 1:length(TXT2)
    t = strsplit(TXT2{i},'T') ;
    b{i} = t{2} ;
end

timeStamps = b;

%still unsure of how to get this to work, the function below doesn't handle
%it properly, see https://groups.google.com/forum/#!topic/bonsai-users/WD6mV94KAQs