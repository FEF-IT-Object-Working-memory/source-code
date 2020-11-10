%% In the name of Allah

%% SVM analysis for cross validating the spike sorting quality 
% Rezayat. et al,  Nat. Comm. 2020
% "Frontotemporal Coordination Predicts Working Memory Performance and its Local Neural Signatures"
% Fig S1f
clear all
close all
clc

%% Load data
load([cd '\Data\FEF_SVM.mat'])
load([cd '\Data\IT_SVM.mat'])

%% Extract Pairs Performance 
val_u=nan;
val_session_fef=[];
for ss=1:size(t_fefT,1)
   val_h=[];
    for uni=1:size(t_fefT,2)
        
        if ~isempty (t_fefT(ss, uni).pt )
          val_h(uni)=(nanmean(t_fefT(ss, uni).pt));
        else
            if uni==1
                val_h(uni)=val_u;
            end
                
        end
        
    end
            
     val_session_fef =[val_session_fef (val_h)]; 

end

ind_pair_fef=~isnan(val_session_fef);
val_session_it=[];
for ss=1:size(t_itT,1)
    
      val_h=[];
    for uni=1:size(t_itT,2)
        if ~isempty (t_itT(ss, uni).pt )
          val_h(uni)=(nanmean(t_itT(ss, uni).pt));
        else
            if uni==1
                val_h(uni)=val_u;
            end
        
        end
        
    end            
    val_session_it =[val_session_it (val_h)]; 

end
ind_pair_it=~isnan(val_session_it);

val_session_fef=val_session_fef(ind_pair_fef);
val_session_it=val_session_it(ind_pair_it);

val_session_fef_time=[];
for ss=1:size(t_fef,1)
   val_h=[];
    for uni=1:size(t_fef,2)
         for ti=1:size(t_fef,3)
        if ~isempty (t_fef(ss, uni,ti).pt )
          val_h(uni,ti)=(nanmean(t_fef(ss, uni,ti).pt));
        else
            if uni==1
                val_h(uni,ti)=val_u;
            end
                
        end
         end
        
    end
            
        val_session_fef_time = [val_session_fef_time;(val_h)]; 
 
end
val_session_it_time=[];
for ss=1:size(t_it,1)
    
      val_h=[];
    for uni=1:size(t_it,2)
       for ti=1:size(t_it,3)
           if ~isempty (t_it(ss, uni,ti).pt )
           
          val_h(uni,ti)=(nanmean(t_it(ss, uni,ti).pt,1));
      
        else
            if uni==1
                val_h(uni,ti)=val_u;
            end
            end
        end
        
    end
          
      val_session_it_time = [val_session_it_time;(val_h)]; 
end


val_session_it_time=val_session_it_time(ind_pair_it,:);
val_session_fef_time=val_session_fef_time(ind_pair_fef,:);

%% Fig. S2a SVM average classifier performance over time

figure('name','Fig. S2a, Rezayat.E, et al')
ax=subplot(1,1,1);
hold on
niceplot(val_session_it_time*100,[1:49],1,1,0.5,0)

niceplot(val_session_fef_time*100,[1:49],1,1,0,1)
ylabel('Classifier performance(%)')
xlabel('Time across recording session (a.u.)')
set(gca,'fontsize',14,'fontweight','bold')
ylim([90 100])
xlim([1 49])
ax.XTick=[1 49];
ax.YTick=[90 100];
text(35,96.5,'IT','color',[1,0.5,0],'fontsize',20)
text(15,96.5,'FEF','color',[1,0,1],'fontsize',20)
export_file_to_excel([[1:49]',nanmean(val_session_fef_time*100)', nanmean(val_session_it_time*100)' ],'figS2a',14)

%% Fig. S2a inset SVM average classifier performance 

figure('name','Fig. S2a Inset, Rezayat.E, et al')
 ax=subplot(121);
[counts,centers]=hist(val_session_fef*100,[90:1:100]);
counts=counts/sum(counts)*100;
bar(centers,counts)
xlabel('Classifier performance (%)')
set(gca,'fontsize',14,'fontweight','bold')
ylabel('% of pairs ')
title ('FEF')
xlim([90 100])
ax.XTick=[90 100];
ax.YTick=[0 80];
ylim([0 80])
export_file_to_excel([counts',centers' ],'figS2afefinssset',9)

ax=subplot(122);
[counts,centers]=hist(val_session_it*100,[90:1:100]);
counts=counts/sum(counts)*100;
bar(centers,counts)
xlabel('Classifier performance (%)')
xlim([90 100])
set(gca,'fontsize',14,'fontweight','bold')
ylabel('% of pairs ')
title ('IT')
ax.XTick=[90 100];
ax.YTick=[0 80];
ylim([0 80])
export_file_to_excel([counts',centers' ],'figS2aItinset',9)

%% Stat for SVM average classifier performance

disp('%%% State for SVM average classifier performance across sessions %%%')
x_h=val_session_fef*100;
p=signrank(x_h,50);
n=length(x_h);
disp (strcat('SVM average classifier performance , FEF; ',' M',' =',num2str(nanmean(x_h)),' + ',...
    num2str(nanstd(x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p),' compared to chance level, i.e. 50%'));

x_h=val_session_it*100;
p=signrank(x_h,50);
n=length(x_h);
disp (strcat('SVM average classifier performance , IT; ',' M',' =',num2str(nanmean(x_h)),' + ',...
    num2str(nanstd(x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p),' compared to chance level, i.e. 50%'));

