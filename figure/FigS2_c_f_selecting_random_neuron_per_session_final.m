%% In the name of Allah

%% Rezayat. et al,  Nat. Comm. 2020
% "Frontotemporal Coordination Predicts Working Memory Performance and its Local Neural Signatures"
% Fig S2 c f
close all
clc
clear all
%% control for one neuron form each site fig 5e
load([cd '\Mat\PPL\control_selection.mat'])

for iteri=1:1000
    idx_control=logical(zeros(length(ind_neurons_resp),1));
    ses=[];
    for ss=1:92
        ind_h= find(ismember(ind_neurons_resp,ss));
        if ~isempty(ind_h)
            idx_control(ind_h(randi(length(ind_h),1)))=1;
        end
    end
    
    x_h_=mi_out_n_(idx_control);
    y_h_=mi_in_n_(idx_control);
    x_h=[];x_h=x_h_(~isnan(x_h_)&~isnan(y_h_));
    y_h=[];y_h=y_h_(~isnan(x_h_)&~isnan(y_h_));
    
    p_mi_in_out(iteri)=signrank(x_h,y_h);
    num_in_out(iteri)=length(y_h);
    mean_val(iteri)=nanmedian(y_h)-nanmedian(x_h);
end


%%
figure('name','Control for Fig. 5e')
bins=0.0:0.01:0.2;
[counts,centers]=hist(mean_val,bins);
counts=counts/sum(counts)*100;
export_file_to_excel([centers',counts'],'figS2c',9)

bar(centers,counts,1)
line([0 0],ylim)
set(gca,'fontsize',14,'fontweight','bold')
 text( 0.05,5,strcat(num2str((sum(p_mi_in_out<0.05)/length(p_mi_in_out)*100)),'% of valuse <0.05'))
signrank(mean_val)
ylabel('%')
xlabel('p values (a.u.)')
 xlim([-0.02 0.2])
title('Control for Fig. 5e (200 times)')
figure('name','Control for Fig. 5e')
hist(num_in_out)
title ('Total 74 from 43 sites')

%%

v_h=mi_in_n_-mi_out_n_;

same_pairs=[];pp=0;
for ss=1:92
    ind_h=[];ind_h= find(ismember(ind_neurons_resp,ss));
    if length(ind_h)>1
        P_h= nchoosek(ind_h,2);
        if size(P_h,1)>0
            if size(P_h,1)<2
                
                ll= size(same_pairs,1);
                same_pairs(ll+1:ll+1,:)= [[v_h(P_h)]];
            else
                
                same_pairs= [same_pairs; v_h(P_h)];
                
            end
        end
    end
end


figure
hold on
x_h_=same_pairs(:,1);
y_h_=same_pairs(:,2);
x_h=[];x_h=x_h_(~isnan(x_h_)&~isnan(y_h_));
y_h=[];y_h=y_h_(~isnan(x_h_)&~isnan(y_h_));
scatter(x_h,y_h)

[r p]=correlation(x_h,y_h,10);
num_in_out(iteri)=length(y_h);
mean_val(iteri)=nanmedian(y_h)-nanmedian(x_h);
x=x_h;y=y_h;
tbl = table(x,y);
lm = fitlm(tbl,'linear');
b=lm.Coefficients{1,1};
m=lm.Coefficients{2,1};
xx=-0.5:0.1:.5;
% [r,p]=corr(x_h,y_h,'type','kendal');
% [r,p]=corr(x_h,y_h,'type','kendal')
cohen_r2=(r^2)

length(x_h)

plot(xx,m*xx+b,'r')


xlabel(' SPl Cr-Wr[neuron1](a.u.)')
ylabel(' SPl Cr-Wr[neuron2](a.u.)')
ylim([-0.6 0.6])
xlim([-0.6 0.6])
% line([0.5 0.5],ylim)
% line(xlim,[0 0])
set(gca,'fontsize',14,'fontweight','bold')
text( 0.10,0.52,strcat('r = ',num2str(r,2)),'fontsize',14,'fontweight','bold')
text( 0.10,0.40,strcat('p = ',num2str(p,2)),'fontsize',14,'fontweight','bold')
text( 0.10,0.28,strcat('n = ',num2str(length(x_h))),'fontsize',14,'fontweight','bold')

%%
ses_index=1:length(ind_neurons_resp);
P_h= nchoosek(ses_index,2);

for iteri=1:100
    ind_h=randperm(size(P_h,1),51);
    val_h= v_h(P_h(ind_h,:));
    x_h=val_h(:,1); y_h=val_h(:,2);
    x_h_=x_h(~(isnan(y_h)&isnan(x_h)));
    y_h_=y_h(~(isnan(y_h)&isnan(x_h)));
    x_h=x_h_;y_h=y_h_;
    [r_a(iteri),p]=correlation(x_h,y_h,10);
    
end
%%
figure
hold on
[rr cc]=hist(r_a,-1:0.1:1);
export_file_to_excel([cc',[rr/sum(rr)]'],'figS2f',9)

bar(cc,rr/sum(rr),'barwidth',1)
line([r r],ylim,'color','r','linewidth',2)
set(gca,'fontsize',14,'fontweight','bold')
xlim([-0.5 0.5])
xlabel('Correlation value (a.u.)')
ylabel('% pairs')
text( r,0.3,strcat('r = ',num2str(r,2)),'fontsize',14,'fontweight','bold')
text( r,0.25,strcat('p = ',num2str(sum(r_a>r)/length(r_a),2)),'fontsize',14,'fontweight','bold')

