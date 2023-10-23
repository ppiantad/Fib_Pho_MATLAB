
%collect all animal names and blockpaths from user, label blockpaths with
%animal name
% animalNames= {
%    'D1-OP-1';...
%    'D1-OP-2';...
%    'D1-iOP-3';... %10 mW
%    'D1-iOP-4';... %10 mW
%    };
%   
% blockpaths= {
%     'g:\MATLAB\TDTbin2mat\Photometry\127 & D1-OP-1\D1-OP-1\D1-OP-1-200227-132526'...
%     'g:\MATLAB\TDTbin2mat\Photometry\128 & D1-OP-2\D1-OP-2\D1-OP-2-200227-134925'...
%     'G:\MATLAB\TDTbin2mat\Photometry\eOP & iOP 10 mW Laser Test (Sensor D)\D1-iOP-3-210623-190822'...
%     'G:\MATLAB\TDTbin2mat\Photometry\eOP & iOP 10 mW Laser Test (Sensor D)\D1-iOP-4-210623-192535'...
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
% 
% 
% %make table with all animal info
% 
% Table=table(animalNames,blockpaths',groups,dataStruct_names,startTime',Channel_405_name,Channel_465_name,whichStreams, Pu1_epoc);
% Headers={'animalNames','blockpath','group','dataStruct_names','startTime','Channel_405_name','Channel_465_name','whichStreams', 'Pu1_epoc'};
% Table.Properties.VariableNames([1:9])=Headers;
% %fill in animal names

%%
%collect all animal names and blockpaths from user, label blockpaths with
%animal name
animalNames= {
   'D2-OP-1';...
   'D2-OP-2';...
   };
  
blockpaths= {
    'E:\MATLAB\TDTbin2mat\Photometry\125 & D2-OP-1\D2-OP-1\D2-OP-1-200227-142348'...
    'E:\MATLAB\TDTbin2mat\Photometry\126 & D2-OP-2\D2-OP-2\D2-OP-2-200227-144752'...
    };


groupNames={'eNpHR','mCherry'};

groups= {'eNpHR';'eNpHR'};
whichStreams=[12;12];
numAnimals=numel(animalNames);
dataStruct_names=animalNames;
startTime=zeros(1,numAnimals);
Channel_405_name=animalNames;
Channel_465_name=animalNames;


%make table with all animal info

Table=table(animalNames,blockpaths',groups,dataStruct_names,startTime',Channel_405_name,Channel_465_name,whichStreams);
Headers={'animalNames','blockpath','group','dataStruct_names','startTime','Channel_405_name','Channel_465_name','whichStreams'};
Table.Properties.VariableNames([1:8])=Headers;
% %fill in animal names

%%

%run TDTbin2mat on all of these, fill in the table, and store the datastructures in workspace
%for each animal....
zall_all={};ts1_all={};tCS_all={};zCS_all={};  
for aaa=1:numAnimals
    animalName=Table.animalNames{aaa};
    blockpath=Table.blockpath{aaa};
    whichStreams=Table.whichStreams(aaa);
    Pu1_epoc=Table.Pu1_epoc(aaa);
   [ts1, zall] = OFLcond_BLto5s_FN(animalName, blockpath, whichStreams,Pu1_epoc);
    ts1_all{aaa,1}=ts1;
    zall_all{aaa,1}=zall;
 
    
       
end
%%
% % ts1_all_mat=cell2mat(ts1_all);
% [rows,~]=size(zall_all);
% cc=[];
% for p=1:rows
%     [~,cc(p)]=size(zall_all{p,1});
% 
% 
% end
% shortest=min(cc);
% zall_all_mat=[];
% for p=1:rows
%     zall_all_mat(p,:)=zall_all{p,1}(1,1:shortest);
% end


%%
zall_mat=cell2mat(zall_all);
clear zall_all %clear zall_all to save space
% tCS_all_mat=cell2mat(tCS_all);



%DATA ORGANIZATION - assuming inputted one group at a time
zall_mean=mean(zall_mat,1);
sd_zall=std(zall_mat,1);
SEM_all = std(zall_mat)/sqrt(size(zall_mat,1));


%%
%if  all matrices are made, clear the cells



%%
figure (1)
hold on
a=0;

    p1=patch([a a+2 a+2 a], [-4 -4 5 5], [0.9100 0.4100 0.1700], 'EdgeColor','none');
    p1.FaceAlpha = .25;
%     p4=patch([a+28 a+30 a+30 a+28], [-4 -4 22 22], [1 .8 .8], 'EdgeColor','none');



lo=zall_mean-SEM_all;
hi=zall_mean+SEM_all;
xxx1=[ts1_all{1,:}, ts1_all{1,:}(end:-1:1)];
yyy1=[lo, hi(end:-1:1)];
hp1= fill(xxx1,yyy1,[ .8 .8 .8]); hold on;
set(hp1,'FaceColor', [ .8 .8 .8],'EdgeColor','none');hold on;
set(gca,'Layer','Top')
plot(ts1_all{1,:},zall_mean);

%%
% lo2=zCS_mean-sd_zCS;
% hi2=zCS_mean+sd_zCS;
% xxx2=[ts1_all{1,:}, ts1_all{1,:}(end:-1:1)];
% yyy2=[lo2, hi2(end:-1:1)];
% 
% figure (2)
% hold on
% 
% 
% p2=patch([0 28 28 0], [-4 -4 22 22], [.8 1 1], 'EdgeColor','none');
% p3=patch([28 30 30 28], [-4 -4 22 22], [1 .8 .8], 'EdgeColor','none');
% 
% 
% 
% 
% hp2= fill(xxx2,yyy2,[ .8 .8 .8]); hold on;
% set(hp2,'FaceColor', [ .8 .8 .8],'EdgeColor','none');hold on;
% set(gca,'Layer','Top')
% 
% plot(tCS_all{1,1}(1,:), zCS_mean);
%%
%SPLIT CSs into bins of 5
z_bin5={};first5={};last5={};z_bin_nBL={};
for bin = 1:4
    %set bin
    a=5*(bin-1)+1;
    z_this_bin=[];z_this_bin_nBL=[];
    %locate animal and compile first [bin size] for each animal
    for mouse=1:numAnimals
        start=20*(mouse-1);
        start=start+a; fin=start+2;
        z_this_bin= [z_this_bin; zall_mat(start:fin,:)];
       
        
        
        
        
        
        
        if bin==1
            zCS1(mouse,:)=zall_mat(start,:);
            first5{mouse,1}=mean(zall_mat(start:fin,:),1);
            peak1_mouse(mouse)=max(first5{mouse,1}(:,ts1(1,:) < 2 & ts1(1,:) > 0));
            zs_bin1=z_this_bin;
            
            
            
%             first5US{mouse,1}=mean(zBL4US_mat(start:fin,:),1);
%             peakUS_mouse(mouse)=max(first5US{mouse,1}(:,ts1(1,:) < 20 & ts1(1,:) > 16));
            
            
            
            %find top 20% activity in each trial
            z_0to2=zall_mat(start:fin,ts1(1,:) < 2 & ts1(1,:) > 0);
            len_z=numel(z_0to2(1,:));
            top20=round(0.2*len_z);
            
%              z_16to20=zBL4US_mat(start:fin,ts1(1,:) < 20 & ts1(1,:) > 16);
%             len_zUS=numel(z_16to20(1,:));
%             top20US=round(0.2*len_zUS);
            
            for trial = 1:5
                sort_z=sort(z_0to2(trial,:),'descend');
                cap(mouse,trial)=mean(sort_z(1:top20));
                
%                 sort_zUS=sort(z_16to20(trial,:),'descend');
%                 capUS(mouse,trial)=mean(sort_zUS(1:top20));
                
            end
            
            
            
        elseif bin ==4
            last5{mouse,1}=mean(zall_mat(start:fin,:),1);
            zs_bin4=z_this_bin;
            
            
            
            peak10CS_mouse(mouse)=max(last5{mouse,1}(:,ts1(1,:) < 2 & ts1(1,:) > 0));
            
%             last5US{mouse,1}=mean(zBL4US_mat(start:fin,:),1);
%             peak10US_mouse(mouse)=max(last5US{mouse,1}(:,ts1(1,:) < 20 & ts1(1,:) > 16));
            
            
            
            
            
             %find top 20% activity in each trial
             z_0to2=zall_mat(start:fin,ts1(1,:) < 2 & ts1(1,:) > 0);
             len_z=numel(z_0to2(1,:));
             top20=round(0.2*len_z);
%             
%             z_16to20=zBL4US_mat(start:fin,ts1(1,:) < 20 & ts1(1,:) > 16);
%             len_zUS=numel(z_16to20(1,:));
%             top20US=round(0.2*len_zUS);
            
            for trial = 1:4
                sort_z=sort(z_0to2(trial,:),'descend');
                cap(mouse,trial+4)=mean(sort_z(1:top20));
                
%                 sort_zUS=sort(z_16to20(trial,:),'descend');
%                 capUS(mouse,trial+4)=mean(sort_zUS(1:top20));
            end
            
            
        end
        
    end
    
    
    %average across animals
    z_bin{bin,1}=mean(z_this_bin,1);

    bin_sd{bin,1}=(std(z_this_bin,1))./((numel(z_this_bin(:,1)))^.5);% sd is actually SEM as of 10/24/2019

    
end
cap_av=[mean(cap(:,1:3),2) mean(cap(:,4:6),2) ];


%%
%Plot binned data

figure (2)
hold on;

for m=1:6
    sp(m)=subplot(6,1,m);
    
    
    
    p5=patch([0 28 28 0], [-4 -4 45 45], [.8 1 1], 'EdgeColor','none');
    p6=patch([28 30 30 28], [-4 -4 45 45], [1 .8 .8], 'EdgeColor','none');
    hold on;
    
    
    
    lo=z_bin{m,1}-bin_sd{m,1};
    hi=z_bin{m,1}+bin_sd{m,1};
    xxx1=[ts1_all{1,:}, ts1_all{1,:}(end:-1:1)];
    yyy1=[lo, hi(end:-1:1)];
    hp1= fill(xxx1,yyy1,[ .8 .8 .8]); hold on;
    set(hp1,'FaceColor', [ .8 .8 .8],'EdgeColor','none');hold on;
    
    plot(ts1,z_bin{m,1}); hold on; 
    t=title(sprintf('CS %d to %d ', (m-1)*5+1, (m-1)*5+6));
    set(t,'position', [-2.1,3.4,0])
    set(gca,'Layer','top')
    if m~=6
        set(gca,'xtick',[])
    end
    ylim([-10 60])

end


%%
%CS 1-3 and 28-30 on same graph

figure (6)
hold on

p5=patch([0 28 28 0], [-4 -4 45 45], [.8 1 1], 'EdgeColor','none');
p6=patch([28 30 30 28], [-4 -4 45 45], [1 .8 .8], 'EdgeColor','none');
hold on


%SD bin 1
lo=z_bin{1,1}-bin_sd{1,1};
hi=z_bin{1,1}+bin_sd{1,1};
xxx1=[ts1_all{1,:}, ts1_all{1,:}(end:-1:1)];
yyy1=[lo, hi(end:-1:1)];
hp1= fill(xxx1,yyy1,[ 0 .5 .1],'facealpha',.5); hold on;
set(hp1,'FaceColor', [0 .5 .1],'EdgeColor','none');hold on;



%SD bin 10
lo=z_bin{10,1}-bin_sd{10,1};
hi=z_bin{10,1}+bin_sd{10,1};
xxx1=[ts1_all{1,:}, ts1_all{1,:}(end:-1:1)];
yyy1=[lo, hi(end:-1:1)];
hp1= fill(xxx1,yyy1,[ .8 .8 .8],'facealpha',.7); hold on;
set(hp1,'FaceColor', [ .8 .8 .8],'EdgeColor','none');hold on;



p1to5=plot(ts1,z_bin{1,1},'color',[0 .5 .1],'LineWidth',2); hold on;


p26to30=plot(ts1,z_bin{10,1},'color',[0 0 0],'LineWidth',2); hold on;
ylim([-10 60])

legend([p1to5 p26to30],'CS 1 to 3','CS 28 to 30','AutoUpdate','off')

set(gca,'Layer','top')
%%
%CS 1-3 and 28-30 on same graph
%exclude second to last trial for 1-3 bin (CS2 for YFP 10)
figure (7)
hold on

p5=patch([0 28 28 0], [-4 -4 45 45], [.8 1 1], 'EdgeColor','none');
p6=patch([28 30 30 28], [-4 -4 45 45], [1 .8 .8], 'EdgeColor','none');
hold on


%modify bin 1
z_bin1_ex=[zs_bin1(1:13,:); zs_bin1(15,:)];
mean_zbin1_ex=mean(z_bin1_ex,1);
bin1_ex_sd=(std(z_bin1_ex,1))./((14)^.5);


z_bin1_ex_nBL=[zs_bin1_nBL(1:13,:); zs_bin1_nBL(15,:)];
mean_zbin1_ex_nBL=mean(z_bin1_ex_nBL,1);
bin1_ex_sd_nBL=(std(z_bin1_ex_nBL,1))./((14)^.5);

%SD bin 1
lo=mean_zbin1_ex-bin1_ex_sd;
hi=mean_zbin1_ex+bin1_ex_sd;
xxx1=[ts1_all{1,:}, ts1_all{1,:}(end:-1:1)];
yyy1=[lo, hi(end:-1:1)];
hp1= fill(xxx1,yyy1,[ 0 .5 .1],'facealpha',.5); hold on;
set(hp1,'FaceColor', [0 .5 .1],'EdgeColor','none');hold on;



%SD bin 10
lo=z_bin{10,1}-bin_sEM_nbl{10,1};
hi=z_bin{10,1}+bin_sEM_nbl{10,1};
xxx1=[ts1_all{1,:}, ts1_all{1,:}(end:-1:1)];
yyy1=[lo, hi(end:-1:1)];
hp1= fill(xxx1,yyy1,[ .8 .8 .8],'facealpha',.7); hold on;
set(hp1,'FaceColor', [ .8 .8 .8],'EdgeColor','none');hold on;



p1to5=plot(ts1,mean_zbin1_ex,'color',[0 .5 .1],'LineWidth',2); hold on;


p26to30=plot(ts1,z_bin{10,1},'color',[0 0 0],'LineWidth',2); hold on;
ylim([-10 60])

legend([p1to5 p26to30],'CS 1 to 3','CS 28 to 30','AutoUpdate','off')

set(gca,'Layer','top')

%%

%Heatmap each animal
figure (4)
start=[];fin=[];

for an = 1:numAnimals
    start=30*(an-1)+1;fin=start+29;
    subplot(6,2,an)
    imagesc(ts1, 1, zall_mat(start:fin,:));
    colormap('jet'); % c1 = colorbar;
    title(sprintf(animalNames{an,1}));
 
    
end
%%
%Color traces each animal
figure (5)
start=[];fin=[];
col=jet(30);
for an = 1:numAnimals
    start=30*(an-1);
    subplot(6,2,an)
    for CS=1:30
        bb=start+CS;
        plot(ts1, zall_mat(bb,:),'Color',col(CS,:));hold on;
    end
    title(sprintf(animalNames{an,1}));
 
    
end

%%
%short traces
figure (9)
hold on
%plot([0;0],[-20; 20],'--','color',[.8 .8 .8],'LineWidth',2);

%SUBPLOT 1: CS

%SD bin 1
% % lo=z_bin{1,1}-bin_sd{1,1};
% % hi=z_bin{1,1}+bin_sd{1,1};
% lo=mean_zbin1_ex-bin1_ex_sd;
% hi=mean_zbin1_ex+bin1_ex_sd;
% xxx1=[ts1_all{1,:}, ts1_all{1,:}(end:-1:1)];
% yyy1=[lo, hi(end:-1:1)];
% hp1= fill(xxx1,yyy1,[.96 .47 .49],'facealpha',.5); hold on; %old green was 0 .5 .1
% set(hp1,'FaceColor', [.96 .47 .49],'EdgeColor','none');hold on;



%SD bin 10
lo=z_bin{10,1}-bin_sd{10,1};
hi=z_bin{10,1}+bin_sd{10,1};
xxx1=[ts1_all{1,:}, ts1_all{1,:}(end:-1:1)];
yyy1=[lo, hi(end:-1:1)];
hp1= fill(xxx1,yyy1,[.56 .02 .03],'facealpha',.5); hold on;
set(hp1,'FaceColor', [.56 .02 .03],'EdgeColor','none');hold on;



% p1to3=plot(ts1,mean_zbin1_ex,'color',[.96 .47 .49],'LineWidth',2);
% p1to3=plot(ts1,z_bin{1,1},'color',[.96 .47 .49],'LineWidth',2); hold on; %light red
p28to30=plot(ts1,z_bin{10,1},'color',[.56 .02 .03],'LineWidth',2); hold on; %dark red

%add scale bar
plot([-8; -7], [0;0],'-k', [0; 0], [7;8],'k', 'LineWidth',2);

% text(-4.2, -1.5, 'z=1', 'HorizontalAlignment','right');
% text(-3.5, -2.5, '1 s', 'HorizontalAlignment','center');
set(gca, 'XTick',[],'XTickLabel',[],'xcolor',[1 1 1]);

xlim([-10 10])
ylim([-10 40])
%%
%SUBPLOT 2:US
figure (10)
hold on;
% plot([28;28],[-20; 20],'--','color',[.8 .8 .8],'LineWidth',2);
% plot([30;30],[-20; 20],'--','color',[.8 .8 .8],'LineWidth',2);


% %SD bin 1
% lo=z_bin_nBL{1,1}-bin_sEM_nbl{1,1};
% hi=z_bin_nBL{1,1}+bin_sEM_nbl{1,1};
% xxx1=[ts1_all{1,:}, ts1_all{1,:}(end:-1:1)];
% yyy1=[lo, hi(end:-1:1)];
% hp1= fill(xxx1,yyy1,[ .95 .63 .2],'facealpha',.5); hold on;
% set(hp1,'FaceColor', [.95 .63 .2],'EdgeColor','none');hold on;



%SD bin 10
lo=z_bin_nBL{10,1}-bin_sEM_nbl{10,1};
hi=z_bin_nBL{10,1}+bin_sEM_nbl{10,1};
xxx1=[ts1_all{1,:}, ts1_all{1,:}(end:-1:1)];
yyy1=[lo, hi(end:-1:1)];
hp1= fill(xxx1,yyy1,[ .98 .51 .05],'facealpha',.7); hold on;
set(hp1,'FaceColor', [ .98 .51 .05],'EdgeColor','none');hold on;

 %p1to3=plot(ts1,z_bin_nBL{1,1},'color',[.95 .63 .2],'LineWidth',2); hold on; %light orange
p28to30=plot(ts1,z_bin_nBL{10,1},'color',[.98 .51 .05],'LineWidth',2); hold on; %dark orange

% plot([27; 27], [-2;-1],'-k', [27;28], [-2;-2],'k', 'LineWidth',2);
% text(26.8, -1.5, 'z=1', 'HorizontalAlignment','right');
% text(27.5, -2.5, '1 s', 'HorizontalAlignment','center');
set(gca, 'XTick',[],'YTick',[],'XTickLabel',[],'YTickLabel',[],'xcolor',[1 1 1],'ycolor',[1 1 1]);


plot([21;22], [0;0],'-k', [30; 30], [8;10],'k', 'LineWidth',2);

xlim([21 37])
ylim([-40 10])%-9 to 11 for yfp
%%
%Exemplar CS
figure (11)
hold on
plot([0;0],[-20; 20],'--','color',[.8 .8 .8],'LineWidth',2);
plot(ts1,zall_mat(89,:),'color',[.96 0 0],'LineWidth',2); hold on;

xlim([-5 5])
ylim([-5 40])%-5 to 10 for YFP



%Exemplar US
figure (12)
hold on;
plot([28;28],[-20; 20],'--','color',[.8 .8 .8],'LineWidth',2);
plot([30;30],[-20; 20],'--','color',[.8 .8 .8],'LineWidth',2);
plot(ts1,zBL4US_mat(77,:),'color',[1 .64 0],'LineWidth',2); hold on;


xlim([23 32])
ylim([-40 5])%-4 15 for YFp



%%

%Exemplar heatmap CS
figure (13); hold on;
imagesc(ts1(ts1>-5 &ts1<5), 1, zall_mat(91:120,(ts1>-5 &ts1<5)));
colormap('jet');  c1 = colorbar;
set(gca,'Visible','off')



%Exemplar heatmap US
figure (14); hold on;
imagesc(ts1(ts1>23 &ts1<32), 1, zall_mat(61:90,(ts1>23 &ts1<32)));
colormap('jet');  c1 = colorbar;
set(gca,'Visible','off')
%%


%all CSs


%average all zs for each animal

start=1;fin=30;
for mouse = 1:numAnimals
    start=30*(mouse-1)+1;fin=start+29;
    zCS_mouse=zall_mat(start:fin,:);
    zCS(mouse,:)=mean(zCS_mouse,1);
    peak_allCS(mouse,1)=max(zCS(mouse,ts1<2 & ts1>0));
  
    
    
    zUS_mouse=zBL4US_mat(start:fin,:);
    zUS(mouse,:)=mean(zUS_mouse,1);
    peak_allUS(mouse,1)=max(zCS(mouse,ts1<30 & ts1>28));
    
    
end
zCS_avg=mean(zCS,1);
zCS_sem=(std(zCS_mouse,1))./((numel(zCS_mouse(:,1)))^.5);



zUS_avg=mean(zUS,1);
zUS_sem=(std(zUS_mouse,1))./((numel(zUS_mouse(:,1)))^.5);








%%

figure (9)

%all CSs
lo=zCS_avg-zCS_sem;
hi=zCS_avg+zCS_sem;
xxx1=[ts1_all{1,:}, ts1_all{1,:}(end:-1:1)];
yyy1=[lo, hi(end:-1:1)];
hp1= fill(xxx1,yyy1,[ 0 .5 .1],'facealpha',.5); hold on;
set(hp1,'FaceColor', [0 .5 .1],'EdgeColor','none');hold on;
plot(ts1,zCS_avg,'color',[0 .5 .1],'LineWidth',2); hold on;
xlim([-5 4])
ylim([-4 6])

figure (10)
%all USs
lo=zUS_avg-zUS_sem;
hi=zUS_avg+zUS_sem;
xxx1=[ts1_all{1,:}, ts1_all{1,:}(end:-1:1)];
yyy1=[lo, hi(end:-1:1)];
hp1= fill(xxx1,yyy1,[ .8 .8 .8],'facealpha',.7); hold on;
set(hp1,'FaceColor', [ .8 .8 .8],'EdgeColor','none');hold on;
plot(ts1,zUS_avg,'color',[0 0 0],'LineWidth',2); hold on;
xlim([23 32])
ylim([-4 6])











%%
%AUCS and peak -  normalized 0 to 5s 

%-2 to 0; 0 to 2; 28 to 30; 30 to 32;

per_on=[-5 0 26 28 30];
per_off=[0 5 28 30 32];
AUC=[];

%zall_mean; bin 1; bin 10
for num=1:numel(per_on)
    AUC(num)=trapz(zall_mean(:,ts1(1,:) < per_off(num) & ts1(1,:) > per_on(num)));
    AUC_bin1(num)=trapz(z_bin{1,1}(:,ts1(1,:) < per_off(num) & ts1(1,:) > per_on(num)));
    AUC_bin10(num)=trapz(z_bin{10,1}(:,ts1(1,:) < per_off(num) & ts1(1,:) > per_on(num)));
  
    
     for mouse = 1:numAnimals
        AUC_first3(mouse,num)=trapz(first3{mouse,1}(:,ts1(1,:) < per_off(num) & ts1(1,:) > per_on(num)));
        AUC_first3US(mouse,num)=trapz(first3US{mouse,1}(:,ts1(1,:) < per_off(num) & ts1(1,:) > per_on(num)));
        AUC_last3(mouse,num)=trapz(last3{mouse,1}(:,ts1(1,:) < per_off(num) & ts1(1,:) > per_on(num)));
        AUC_last3US(mouse,num)=trapz(last3US{mouse,1}(:,ts1(1,:) < per_off(num) & ts1(1,:) > per_on(num)));
        
        
        %AUC_allCS(mouse,num)=trapz(zCS(mouse,ts1(1,:) < per_off(num) & ts1(1,:) > per_on(num)));
        %AUC_allUS(mouse,num)=trapz(zUS(mouse,ts1(1,:) < per_off(num) & ts1(1,:) > per_on(num)));
                
    end
    
end


%data for excluded trial for YFP 10
z_Y10=mean([zs_bin1(13,:); zs_bin1(15,:)],1);
z_Y10b10=mean([zs_bin10(13,:); zs_bin10(15,:)],1);
AUC_first3_YFP10=trapz(z_Y10(:,ts1(1,:) < 0 & ts1(1,:) > -5));
AUC_firt3_YFP10_0to5=trapz(z_Y10(:,ts1(1,:) < 5 & ts1(1,:) > 0));


AUC_YFP10_US(1,1)=trapz(mean_zbin1_ex_nBL(:,ts1(1,:) < 28 & ts1(1,:) > 26));
AUC_YFP10_US(1,2)=trapz(mean_zbin1_ex_nBL(:,ts1(1,:) < 30 & ts1(1,:) > 28));
AUC_YFP10_US(1,3)=trapz(mean_zbin1_ex_nBL(:,ts1(1,:) < 32 & ts1(1,:) > 30));

peak1_YFP10=max(z_Y10(:,ts1(1,:) < 2 & ts1(1,:) > 0));
peak10_YFP10=max(z_Y10b10(:,ts1(1,:) < 30 & ts1(1,:) > 28));
%Peak - as defined by max in first 2 seconds of CS
for bin = 1:6
   peak(bin)=max( z_bin{bin,1}(:,ts1(1,:) < 2 & ts1(1,:) > 0));
   
end

%%
%AUCS and peak -  normalized 23 to 28s (nbl = new baseline)

%-2 to 0; 0 to 2; 28 to 30; 30 to 32;

per_on=[26 28 30];
per_off=[28 30 32];
AUC=[];

%zall_mean; bin 1; bin 10
for num=1:numel(per_on)
    for mouse = 1:numAnimals
        AUC_first5nbl(mouse,num)=trapz(first3US{mouse,1}(:,ts1(1,:) < per_off(num) & ts1(1,:) > per_on(num)));
        AUC_last5nbl(mouse,num)=trapz(last3US{mouse,1}(:,ts1(1,:) < per_off(num) & ts1(1,:) > per_on(num)));
    end
    
end



%%
%Cumulative distribution of first five seconds


[u,w]=find(ts1(1,:)>0 & ts1(1,:) < 2);
time_in_2=ts1(w);


for bin = 1:6
    AUC_5=trapz(time_in_5,z_bin{bin,1}(1,w));
    for n=1:numel(time_in_5)-1
        next=5087+n;
        cd_x(bin,n)=trapz(time_in_5(1:n+1) ,(z_bin{bin,1}(1,5087:next)));
        cd_p(bin,n)=cd_x(bin,n)/AUC_5;
        
    end
end
%%
%PLOT Cumulative distribution

figure (4)
hold on;

for m=1:6
    sr(m)=subplot(6,1,m);
    
    hold on;
    plot(time_in_5(1:2033),cd_p(m,:)); hold on;
    set(gca,'Layer','top')
    if m~=6
        set(gca,'xtick',[])
    end
   ylim([0 1.1])
    xlim([0 2])

end


%%
%Just CS1

zCS1_avg=mean(zCS1,1);

figure (3)
hold on
patch([0 28 28 0], [-4 -4 45 45], [.8 1 1], 'EdgeColor','none');
patch([28 30 30 28], [-4 -4 45 45], [1 .8 .8], 'EdgeColor','none');
plot(ts1, zCS1_avg);



%%
% %Make each animals' data same length
% %%
% %find smallest 
% sizes=cellfun(@numel,zAll);
% [minLengthAn,which]=min(sizes);
% [~,minCols]=size(zAll{1,which});
% 
% [~,cellCols]=size(zAll);
% 
% for col=1:cellCols
%     zAll{1,col}=zAll{1,col}(:,1:minCols);
%        
% end




% %% FIGURE 1: group average for whole session
% ArchColor = [1 .5 0]; YFPColor = [0 0 0];
%     
% %change patch depending on what you're graphing!
% 
% figure (1); hold on;
% 
% p1=patch([-4 0 0 -4], [-30 -30 30 30], [.8 1 1], 'EdgeColor','none');
% vals = zAll{1,mouse}(sub,:);
% plot(ts1,vals,'Color',whichColor); hold on;
% 
% 
% 
% set(gca, 'Layer', 'top')
    
