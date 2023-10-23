%% set drive
root = upper('G:');

%% Individual animal test
% animalNames= {
%    'D1-eOP-6';...
%    
%    };
%   
% blockpaths= {
%     strcat(root,'\MATLAB\TDTbin2mat\Photometry\eOP & iOP 10 mW Laser Test (Sensor D)\D1-eOP-6-210623-175117')...
%     
%     };
% 
% 
% groupNames={'ChrimsonR'};
% 
% groups= {'ChrimsonR'};
% whichStreams=[12];
% Pu1_epoc = [5];
% numAnimals=numel(animalNames);
% dataStruct_names=animalNames;
% startTime=zeros(1,numAnimals);
% Channel_405_name=animalNames;
% Channel_465_name=animalNames;




%% D1: 10 mW 2s Test Sensor G
% animalNames= {
%    'D1-eOP-1';...
%    'D1-eOP-2';...
%    'D1-eOP-6';...
%    'D1-eOP-7';...
%    };
%   
% blockpaths= {
%     strcat(root,'\MATLAB\TDTbin2mat\Photometry\Laser Test Sensor G\D1-eOP-1-201014-125552')...
%     strcat(root,'\MATLAB\TDTbin2mat\Photometry\Laser Test Sensor G\D1-eOP-2-201014-131416')...
%     strcat(root,'\MATLAB\TDTbin2mat\Photometry\eOP & iOP 10 mW Laser Test (Sensor D)\D1-eOP-6-210623-175117')...
%     strcat(root,'\MATLAB\TDTbin2mat\Photometry\eOP & iOP 10 mW Laser Test (Sensor D)\D1-eOP-7-210623-181010')...
%     };
% 
% 
% groupNames={'ChrimsonR';'ChrimsonR';'ChrimsonR';'ChrimsonR'};
% 
% groups= {'ChrimsonR';'ChrimsonR';'ChrimsonR';'ChrimsonR'};
% whichStreams=[12;12;12;12];
% Pu1_epoc = [3; 3; 5; 5];
% numAnimals=numel(animalNames);
% dataStruct_names=animalNames;
% startTime=zeros(1,numAnimals);
% Channel_405_name=animalNames;
% Channel_465_name=animalNames;

%% D2: 10 mW 2s Test Sensor G
% animalNames= {
%    'D2-eOP-1';...
%    'D2-eOP-6';...
%    'D2-eOP-7';...
%    };
%   
% blockpaths= {
%     strcat(root,'\MATLAB\TDTbin2mat\Photometry\Laser Test Sensor G\D2-eOP-1-201014-133256')...
%     strcat(root,'\MATLAB\TDTbin2mat\Photometry\eOP & iOP 10 mW Laser Test (Sensor D)\D2-eOP-6-210623-171501')...
%     strcat(root,'\MATLAB\TDTbin2mat\Photometry\eOP & iOP 10 mW Laser Test (Sensor D)\D2-eOP-7-210623-173302')...
%     };
% 
% 
% groupNames={'ChrimsonR';'ChrimsonR';'ChrimsonR'};
% 
% groups= {'ChrimsonR';'ChrimsonR';'ChrimsonR'};
% whichStreams=[12; 12; 12];
% Pu1_epoc = [3; 5; 5];
% numAnimals=numel(animalNames);
% dataStruct_names=animalNames;
% startTime=zeros(1,numAnimals);
% Channel_405_name=animalNames;
% Channel_465_name=animalNames;

%% D2: 10 mW 2s Test Sensor G
% D2 MICE WITH NO SIGNAL (MAYBE?)

% animalNames= {
%    'D2-eOP-2';...
%    'D2-eOP-3';...
%    };
%   
% blockpaths= {
%     strcat(root,'\MATLAB\TDTbin2mat\Photometry\Laser Test Sensor G\D2-eOP-2-201014-135055')...
%     strcat(root,'\MATLAB\TDTbin2mat\Photometry\Laser Test Sensor G\D2-eOP-3-201014-140949')...
%     };
% 
% 
% groupNames={'ChrimsonR';'ChrimsonR'};
% 
% groups= {'ChrimsonR';'ChrimsonR'};
% whichStreams=[12;12];
% Pu1_epoc = [3, 3];
% numAnimals=numel(animalNames);
% dataStruct_names=animalNames;
% startTime=zeros(1,numAnimals);
% Channel_405_name=animalNames;
% Channel_465_name=animalNames;

%% D1: 10 mW 2s Test Sensor G
% D1 MICE WITH NO SIGNAL (MAYBE?)

% animalNames= {
%    'D1-eOP-3';...
%    
%    };
%   
% blockpaths= {
%     'strcat(root,'\MATLAB\TDTbin2mat\Photometry\Laser Test Sensor G\D1-eOP-3-201014-123521')...
%     
%     };
% 
% 
% groupNames={'ChrimsonR'};
% 
% groups= {'ChrimsonR'};
% whichStreams=[12];
% Pu1_epoc = [3];
% numAnimals=numel(animalNames);
% dataStruct_names=animalNames;
% startTime=zeros(1,numAnimals);
% Channel_405_name=animalNames;
% Channel_465_name=animalNames;

%%
%Halorhodopsin D1
% animalNames= {
%    'D1-OP-1';...
%    'D1-OP-2';...
%    'D1-iOP-3';... %10 mW
%    'D1-iOP-4';... %10 mW
%    };
%   
% blockpaths= {
%     strcat(root,'\MATLAB\TDTbin2mat\Photometry\127 & D1-OP-1\D1-OP-1\D1-OP-1-200227-132526')...
%     strcat(root,'\MATLAB\TDTbin2mat\Photometry\128 & D1-OP-2\D1-OP-2\D1-OP-2-200227-134925')...
%     strcat(root,'\MATLAB\TDTbin2mat\Photometry\eOP & iOP 10 mW Laser Test (Sensor D)\D1-iOP-3-210623-190822')...
%     strcat(root,'\MATLAB\TDTbin2mat\Photometry\eOP & iOP 10 mW Laser Test (Sensor D)\D1-iOP-4-210623-192535')...
%     };
% 
% 
% groupNames={'eNpHR','eNpHR','eNpHR','eNpHR'};
% 
% groups= {'eNpHR';'eNpHR';'eNpHR';'eNpHR'};
% whichStreams=[12;12;12;12];
% Pu1_epoc = [3; 3; 5; 5];
% numAnimals=numel(animalNames);
% dataStruct_names=animalNames;
% startTime=zeros(1,numAnimals);
% Channel_405_name=animalNames;
% Channel_465_name=animalNames;

%% Halorhodopsin D2
animalNames= {
   'D2-OP-1';...
   'D2-OP-2';...
   };
  
blockpaths= {
    strcat(root,'\MATLAB\TDTbin2mat\Photometry\125 & D2-OP-1\D2-OP-1\D2-OP-1-200227-142348')...
    strcat(root,'\MATLAB\TDTbin2mat\Photometry\126 & D2-OP-2\D2-OP-2\D2-OP-2-200227-144752')...
    };



groupNames={'eNpHR','mCherry'};

groups= {'eNpHR';'eNpHR'};
whichStreams=[12;12];
Pu1_epoc = [3; 3];
numAnimals=numel(animalNames);
dataStruct_names=animalNames;
startTime=zeros(1,numAnimals);
Channel_405_name=animalNames;
Channel_465_name=animalNames;


%% D1: 4 mW 2s Test Sensor G

% animalNames= {
%    'D1-eOP-1';...
%    'D1-eOP-2';...
%    };
%   
% blockpaths= {
%     strcat(root,'\MATLAB\TDTbin2mat\Photometry\4 mW Laser Test Sensor G\D1-eOP-1-201020-115325')...
%     strcat(root,'\MATLAB\TDTbin2mat\Photometry\4 mW Laser Test Sensor G\D1-eOP-2-201020-121732')...
%     };
% 
% 
% groupNames={'ChrimsonR';'ChrimsonR'};
% 
% groups= {'ChrimsonR';'ChrimsonR'};
% whichStreams=[12;12];
% Pu1_epoc = [3, 3];
% numAnimals=numel(animalNames);
% dataStruct_names=animalNames;
% startTime=zeros(1,numAnimals);
% Channel_405_name=animalNames;
% Channel_465_name=animalNames;

%% D2: 4 mW 2s Test Sensor G

% animalNames= {
%    'D2-eOP-1';...
%   
%    };
%   
% blockpaths= {
%     strcat(root,'\MATLAB\TDTbin2mat\Photometry\4 mW Laser Test Sensor G\D2-eOP-1-201020-120450')...
%    
%     };
% 
% 
% groupNames={'ChrimsonR'};
% 
% groups= {'ChrimsonR'};
% whichStreams=[12];
% Pu1_epoc = [3];
% numAnimals=numel(animalNames);
% dataStruct_names=animalNames;
% startTime=zeros(1,numAnimals);
% Channel_405_name=animalNames;
% Channel_465_name=animalNames;

%% D1: 2 mW 20 hz 30s Test Sensor D
%
% animalNames= {
%    'D1-eOP-1';...
%    'D1-eOP-2';...
%    };
%   
% blockpaths= {
%     strcat(root,'\MATLAB\TDTbin2mat\Photometry\2 mW Laser Test Sensor D\143_D1-eOP-1-201026-150529')...
%     strcat(root,'\MATLAB\TDTbin2mat\Photometry\2 mW Laser Test Sensor D\D1-eOP-2-201026-154202')...
%     };
% 
% 
% groupNames={'ChrimsonR';'ChrimsonR'};
% 
% groups= {'ChrimsonR';'ChrimsonR'};
% whichStreams=[12; 12];
% numAnimals=numel(animalNames);
% dataStruct_names=animalNames;
% startTime=zeros(1,numAnimals);
% Channel_405_name=animalNames;
% Channel_465_name=animalNames;

%% D2: 2 mW 20 Hz 30s Test Sensor D
%
% animalNames= {
%    'D2-eOP-1';...
%    };
%   
% blockpaths= {
%     strcat(root,'\MATLAB\TDTbin2mat\Photometry\2 mW Laser Test Sensor D\D2-eOP-1-201026-162458')...
% 
%     };
% 
% 
% groupNames={'ChrimsonR'};
% 
% groups= {'ChrimsonR'};
% whichStreams=[12];
% numAnimals=numel(animalNames);
% dataStruct_names=animalNames;
% startTime=zeros(1,numAnimals);
% Channel_405_name=animalNames;
% Channel_465_name=animalNames;

%% D1: 2 mW 2s on 2s OFF 20 Hz 30s Sensor D
% animalNames= {
%    'D1-eOP-1';...
%    'D1-eOP-2';...
%    };
%   
% blockpaths= {
%     strcat(root,'\MATLAB\TDTbin2mat\Photometry\LASER_TEST_2MW_2S_ON_2S_OFF\D1-eOP-1-201027-124241')...
%     strcat(root,'\MATLAB\TDTbin2mat\Photometry\LASER_TEST_2MW_2S_ON_2S_OFF\D1-eOP-2-201027-131831')...
%     };
% 
% 
% groupNames={'ChrimsonR';'ChrimsonR'};
% 
% groups= {'ChrimsonR';'ChrimsonR'};
% whichStreams=[12;12];
% numAnimals=numel(animalNames);
% dataStruct_names=animalNames;
% startTime=zeros(1,numAnimals);
% Channel_405_name=animalNames;
% Channel_465_name=animalNames;

%% D2: 2 mW 2s on 2s OFF 20 Hz 30s Sensor D
% animalNames= {
%    'D2-eOP-1';...
% 
%    };
%   
% blockpaths= {
%     strcat(root,'\MATLAB\TDTbin2mat\Photometry\LASER_TEST_2MW_2S_ON_2S_OFF\D2-eOP-1-201027-135622')...
%    
%     };
% 
% 
% groupNames={'ChrimsonR'};
% 
% groups= {'ChrimsonR'};
% whichStreams=[12];
% numAnimals=numel(animalNames);
% dataStruct_names=animalNames;
% startTime=zeros(1,numAnimals);
% Channel_405_name=animalNames;
% Channel_465_name=animalNames;

%% D1: 2 mW constant 20 Hz 30s Sensor D
% animalNames= {
%    'D1-eOP-1';...
%    
%    };
%   
% blockpaths= {
%     strcat(root,'\MATLAB\TDTbin2mat\Photometry\143 & D1-eOP-1\143_D1-eOP-1-201030-103606')...
%     
%     };
% 
% 
% groupNames={'ChrimsonR'};
% 
% groups= {'ChrimsonR'};
% whichStreams=[12];
% numAnimals=numel(animalNames);
% dataStruct_names=animalNames;
% startTime=zeros(1,numAnimals);
% Channel_405_name=animalNames;
% Channel_465_name=animalNames;


%%
%make table with all animal info

Table=table(animalNames,blockpaths',groups,dataStruct_names,startTime',Channel_405_name,Channel_465_name,whichStreams,Pu1_epoc);
Headers={'animalNames','blockpath','group','dataStruct_names','startTime','Channel_405_name','Channel_465_name','whichStreams','Pu1_epoc'};
Table.Properties.VariableNames([1:9])=Headers;
% %fill in animal names

zall_all={};ts1_all={};tCS_all={};zCS_all={};  
for aaa=1:numAnimals
    animalName=Table.animalNames{aaa};
    blockpath=Table.blockpath{aaa};
    whichStreams=Table.whichStreams(aaa);
    Pu1_epoc=Table.Pu1_epoc(aaa);
   [ts1, zall] = ChrimsonR_Photometry_test_fn(animalName, blockpath, whichStreams, Pu1_epoc);
    ts1_all{aaa,1}=ts1;
    zall_all{aaa,1}=zall;
 
    
       
end

%%
zall_mat=cell2mat(zall_all);
clear zall_all %clear zall_all to save space
% tCS_all_mat=cell2mat(tCS_all);



%DATA ORGANIZATION - assuming inputted one group at a time
zall_mean=mean(zall_mat,1);
sd_zall=std(zall_mat,1);
SEM_all = std(zall_mat)/sqrt(size(zall_mat,1));

%%
figure (1)
hold on
a=0;
lo=zall_mean-SEM_all;
hi=zall_mean+SEM_all;
round_hi = round(max(hi));
    p1=patch([a a+2 a+2 a], [-4 -4 round_hi round_hi], [0.9100 0.4100 0.1700], 'EdgeColor','none');
    p1.FaceAlpha = .25;
%     p4=patch([a+28 a+30 a+30 a+28], [-4 -4 22 22], [1 .8 .8], 'EdgeColor','none');




xxx1=[ts1_all{1,:}, ts1_all{1,:}(end:-1:1)];
yyy1=[lo, hi(end:-1:1)];
hp1= fill(xxx1,yyy1,[ .8 .8 .8]); hold on;
set(hp1,'FaceColor', [ .8 .8 .8],'EdgeColor','none');hold on;
set(gca,'Layer','Top')
plot(ts1_all{1,:},zall_mean);

%%
figure(3);
    imagesc(ts1, 1, zall_mat(1:20,:));
    colormap('jet'); % c1 = colorbar;
    
figure(4);
    imagesc(ts1, 1, zall_mat(21:41,:));
    colormap('jet'); % c1 = colorbar;
    
%%
hold off;

first_5 = zall_mat([1:5,21:25],:);
second_5 = zall_mat([6:10,26:30],:);
third_5 = zall_mat([11:15,21:35],:);
fourth_5 = zall_mat([16:20,36:40],:);

first_5_mean = mean(first_5);
second_5_mean = mean(second_5);
third_5_mean = mean(third_5);
fourth_5_mean = mean(fourth_5);

plot(ts1, first_5_mean,'DisplayName','Trials 1-5');hold on;plot(ts1, second_5_mean,'DisplayName','Trials 6-10');plot(ts1, third_5_mean,'DisplayName','Trials 11-15');plot(ts1, fourth_5_mean,'DisplayName','Trials 16-20');hold off;

%%
% hold off;
% first_5 = zall_mat([1:5],:);
% second_5 = zall_mat([6:10],:);
% third_5 = zall_mat([11:15],:);
% fourth_5 = zall_mat([16:20],:);
% 
% first_5_mean = mean(first_5);
% second_5_mean = mean(second_5);
% third_5_mean = mean(third_5);
% fourth_5_mean = mean(fourth_5);
% 
% plot(ts1, first_5_mean,'DisplayName','Trials 1-5');hold on;plot(ts1, second_5_mean,'DisplayName','Trials 6-10');plot(ts1, third_5_mean,'DisplayName','Trials 11-15');plot(ts1, fourth_5_mean,'DisplayName','Trials 16-20');hold off;

%%
hold off;

first_5 = zall_mat([1:5,11:15],:);
second_5 = zall_mat([6:10,16:20],:);


first_5_mean = mean(first_5);
second_5_mean = mean(second_5);


plot(first_5_mean,'DisplayName','Trials 1-5');hold on;plot(second_5_mean,'DisplayName','Trials 6-10');hold off;
