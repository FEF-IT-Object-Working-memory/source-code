function Fig1b_FigS1_Performance(resp)
gen_eye(1).X=resp(1).eye.X;gen_eye(1).Y=resp(1).eye.Y;
gen_eye(1).bhv{1}=resp(1).bhv;gen_eye(1).condition{1}=resp(1).condition';
gen_session_inf(1).bhv{1}=resp(1).bhv;gen_session_inf(1).condition{1}=resp(1).condition';
gen_session_inf(1).preff_obj=resp(1).pref_obj_mu;
gen_session_inf(1).npreff_obj=resp(1).npref_obj_mu;
gen_session_inf(1).FEF_artifact=resp(1).FEF_artifact;
gen_session_inf(1).IT_artifact=resp(1).IT_artifact;

load('F:\Data\Reaserch\Thesis\FEF_IT project\Sample_session\sample_locations.mat')

%% Calculation of Performance
trial_num_conditions = [];performance_of_conditions = [];
for ss =1
    obj_ix_h=[];obj_ix_h = [gen_session_inf(ss).preff_obj,gen_session_inf(ss).npreff_obj];
    obj_ix_h = [obj_ix_h(1) setdiff([1:3],obj_ix_h) obj_ix_h(2)];
    for condi =1:6
        if condi<4; obj_h = obj_ix_h(condi);else obj_h=obj_ix_h(condi-3)+3;end
        ix_h = [];
        
        ix_h = ismember(gen_session_inf(ss).condition{:},obj_h)&(gen_session_inf(ss).bhv{:}==1)&(~gen_session_inf(ss).FEF_artifact)&(~gen_session_inf(ss).IT_artifact);
        
        trial_num_conditions.cr(ss,condi) = sum(ix_h);
        ix_h = [];
        
        ix_h = ismember(gen_session_inf(ss).condition{:},obj_h)&(gen_session_inf(ss).bhv{:}==0)&(~gen_session_inf(ss).FEF_artifact)&(~gen_session_inf(ss).IT_artifact);
        
        trial_num_conditions.wr(ss,condi) = sum(ix_h);
    end
end
performance_of_conditions = (trial_num_conditions.cr)./(trial_num_conditions.cr+trial_num_conditions.wr)*100;

performanceT=nanmean(performance_of_conditions,2);                                        % Average of perfomance across 6 conditions
performance_In=nanmean(performance_of_conditions(:,1:3),2);                                        % Average of perfomance across 6 conditions
performance_Out=nanmean(performance_of_conditions(:,4:6),2);                                        % Average of perfomance across 6 conditions

%% Calculation of saccades parameters
Saccade=[];
sac_info(1).saccades=resp(1).saccades;
for ss=1
    sac_time=[];sac_loc_land=[];sac_loc_start=[];sac_speed=[];sac_r=[];sac_theta=[];
    for tri=1:size(gen_eye(ss).X,1)
        if tri<=size(sac_info(ss).saccades,2)
            if ~isempty(sac_info(ss).saccades{1,tri})
                sac_t_h=[]; sac_t_h=[sac_info(ss).saccades{1,tri}.start_time];
                sac_speed_h=[];sac_speed_h=([sac_info(ss).saccades{1,tri}.r])./([sac_info(ss).saccades{1,tri}.stop_time]- [sac_info(ss).saccades{1,tri}.start_time])*1000;
                sac_r_h=[]; sac_r_h=[sac_info(ss).saccades{1,tri}.r];
                sac_theta_h=[]; sac_theta_h=[sac_info(ss).saccades{1,tri}.theta];
                sac_loc_land_h=[];
                sac_loc_start_h=[];
                for saci=1:length(sac_t_h)
                    sac_loc_land_h(saci,:)=sac_info(ss).saccades{1,tri}(saci).stop_pos;
                    sac_loc_start_h(saci,:)=sac_info(ss).saccades{1,tri}(saci).start_pos;
                end
                
                if  ~isempty(sac_t_h)
                    sac_time_h=sac_t_h(find(sac_t_h>1900&sac_t_h<2500))-1800;
                    sac_loc_h=sac_loc_land_h(find(sac_t_h>1900&sac_t_h<2500),:);
                    sac_loc_st_h=sac_loc_start_h(find(sac_t_h>1900&sac_t_h<2500),:);
                    
                    sac_speed_h=sac_speed_h(find(sac_t_h>1900&sac_t_h<2500));
                    sac_r_h=sac_r_h(find(sac_t_h>1900&sac_t_h<2500));
                    sac_theta_h=sac_theta_h(find(sac_t_h>1900&sac_t_h<2500));
                    
                    if  ~isempty(sac_time_h)
                        
                        sac_time(tri)=sac_time_h(1);
                        sac_loc_land(tri,:)=sac_loc_h(1,:);
                        sac_loc_start(tri,:)=sac_loc_st_h(1,:);
                        sac_speed(tri)=sac_speed_h(1);
                        sac_r(tri)=sac_r_h(1);
                        sac_theta(tri)=sac_theta_h(1);
                    else
                        sac_time(tri)=nan;
                        sac_loc_land(tri,:)=[nan nan];
                        sac_loc_start(tri,:)=[nan nan];
                        sac_speed(tri)=nan;
                        sac_r(tri)=nan;
                        sac_theta(tri)=nan;
                    end
                end
            end
        else
            
            sac_time(tri)=nan;
            sac_loc_land(tri,:)=[nan nan];
            sac_loc_start(tri,:)=[nan nan];
            sac_speed(tri)=nan;
            sac_r(tri)=nan;
            sac_theta(tri)=nan;
            
        end
    end
    Saccade.RT{ss} =sac_time;
    Saccade.sac_land{ss}=sac_loc_land ;
    Saccade.sac_start{ss}=sac_loc_start ;
    
    Saccade.sac_speed{ss}=sac_speed;
    Saccade.sac_r{ss}=sac_r;
    Saccade.sac_theta{ss}=sac_theta;
    
    
end

%% Calculation of Reaction Time
RT_cr = [];RT_wr = [];min_RT=50;
for ss=1
    
    ix_h = [];
    ix_h = (gen_session_inf(ss).bhv{:}==1);
    val_h = []; val_h =Saccade.RT{ss}(ix_h);
    val_h(val_h<min_RT)=nan;
    RT_cr(ss) = nanmedian(val_h);
    
    ix_h = [];
    ix_h = (gen_session_inf(ss).bhv{:}==0);
    val_h = []; val_h =Saccade.RT{ss}(ix_h);
    val_h(val_h<min_RT)=nan;
    RT_wr(ss) = nanmedian(val_h);
    
end
RT_cr(RT_cr<min_RT)=nan;RT_wr(RT_wr<min_RT)=nan;

%% Calculation of Saccade Landing Points Variability
SL_cr = [];SL_wr = [];
for ss=1
    right_side=[];left_side=[];
    right_side=((Saccade.sac_theta{ss}*180/pi)>=0)';
    left_side=((Saccade.sac_theta{ss}*180/pi)<0)';
    if ismember(ss,[64 65 78 83])
        right_side=[];left_side=[];
        right_side=(Saccade.sac_land{ss}(:,1)>=0)';
        left_side=(Saccade.sac_land{ss}(:,2)<0)';
    end
    C_m_r=nanmean(Saccade.sac_land{ss}(right_side,:));
    C_m_l=nanmean(Saccade.sac_land{ss}(left_side,:));
    C_m_r_=C_m_r;C_m_l_=C_m_l;
    
    ix_h = []; ix_h = (gen_session_inf(ss).bhv{:}==1);
    val_right=(sqrt(sum((Saccade.sac_land{ss}(right_side&ix_h,:)-repmat(C_m_r,sum(right_side&ix_h),1)).^2,2)));
    val_right=val_right(val_right<10);
    val_right=std( val_right);
    val_left=(sqrt(sum((Saccade.sac_land{ss}(left_side&ix_h,:)-repmat(C_m_l,sum(left_side&ix_h),1)).^2,2)));
    val_left=val_left(val_left<10);
    val_left=std( val_left);
    SL_cr(ss) =nanmean([val_right val_left]);
    
    ix_h = []; ix_h = (gen_session_inf(ss).bhv{:}==0);
    val_right=(sqrt(sum((Saccade.sac_land{ss}(right_side&ix_h,:)-repmat(C_m_r,sum(right_side&ix_h),1)).^2,2)));
    val_right=val_right(val_right<10);
    val_right=std( val_right);
    val_left=(sqrt(sum((Saccade.sac_land{ss}(left_side&ix_h,:)-repmat(C_m_l,sum(left_side&ix_h),1)).^2,2)));
    val_left=val_left(val_left<10);
    val_left=std( val_left);
    SL_wr(ss) =nanmean([val_right val_left]);
    
end

%% Calculation of Microsaccade and Eye Position Parameters
sac_info(1).microsaccades=resp(1).microsaccades;
win_=20;std_=65;psth=@(x) SmoothN(1000*mean(x,1), win_);
micsac_cr=[];mic_direction_cr=[];mic_length_cr=[];eye_pos_cr=[];
micsac_wr=[];mic_direction_wr=[];mic_length_wr=[];eye_pos_wr=[];
for ss=1
    mic_h=zeros(size(gen_eye(ss).X,1),2100);
    mic_dir_h=zeros(size(gen_eye(ss).X,1),2100);
    mic_r_h=zeros(size(gen_eye(ss).X,1),2100);
    for tri=1:size(sac_info(ss).microsaccades,2)
        if ~isempty(sac_info(ss).microsaccades{1,tri})
            times=[sac_info(ss).microsaccades{1,tri}.time];
            times=times(times<2100);
            
            mic_h(tri,times)=1;
            mic_dir_h(tri,times)=sac_info(ss).microsaccades{1,tri}.theta;
            mic_r_h(tri,times)=sac_info(ss).microsaccades{1,tri}.r;
            
        end
    end
    for c=1:6
        ix_h_cr=ismember(gen_eye(ss).condition{:},c)&(gen_eye(ss).bhv{:}==1);
        ix_h_wr=ismember(gen_eye(ss).condition{:},c)&(gen_eye(ss).bhv{:}==0);
        
        micsac_cr(ss,c,:)=psth(mic_h(ix_h_cr,:));
        micsac_wr(ss,c,:)=psth(mic_h(ix_h_wr,:));
        
        mic_direction_cr(ss,c,:)=circ_mean(mic_dir_h(ix_h_cr,:));
        mic_direction_wr(ss,c,:)=circ_mean(mic_dir_h(ix_h_wr,:));
        
        mic_length_cr(ss,c,:)=nanmean(mic_r_h(ix_h_cr,:),1);
        mic_length_wr(ss,c,:)=nanmean(mic_r_h(ix_h_wr,:),1);
        
        eye_pos_cr(ss,c,:,:)=[nanmean(gen_eye(ss).X(ix_h_cr,1:2100),1);nanmean(gen_eye(ss).Y(ix_h_cr,1:2100),1)];
        eye_pos_wr(ss,c,:,:)=[nanmean(gen_eye(ss).X(ix_h_wr,1:2100),1);nanmean(gen_eye(ss).Y(ix_h_wr,1:2100),1)];
        
    end
    
end

%% Fig. 1b, Mean performance of individual monkeys
ind_ses=1;                                                               % All sessions  92
ind_ses_m1=setdiff(ind_ses,[52:92]);                                        % Sessions from 1-51 is for Monkey 1
ind_ses_m2=setdiff(ind_ses,[1:51]);                                         % Sessions from 52-92 is for Monkey 2

figure('name','Fig. 1b, Rezayat.E, et al')
ax=subplot(1,1,1);
hold on
bar([1],nanmean(performanceT(ind_ses)),0.2,'FaceColor',[0.1 0.1 0.1])
bar([1]+0.2,nanmean(performanceT(ind_ses_m1)),0.2,'FaceColor',[0.5 0.5 0.5])
bar([1]+0.4,nanmean(performanceT(ind_ses_m2)),0.2,'FaceColor',[0.3 0.3 0.3])

val_h=[]; val_h=performanceT(ind_ses);
errorbar([1],nanmean(val_h),nanstd(val_h)/sqrt(length(val_h)),'color',[0.1 0.1 0.1],'linewidth',2)

val_h=[]; val_h=performanceT(ind_ses_m1);
errorbar([1]+0.2,nanmean(val_h),nanstd(val_h)/sqrt(length(val_h)),'color',[0.5 0.5 0.5],'linewidth',2)

val_h=[]; val_h=performanceT(ind_ses_m2);
errorbar([1]+0.4,nanmean(val_h),nanstd(val_h)/sqrt(length(val_h)),'color',[0.3 0.3 0.3],'linewidth',2)

ylim([50 80])
ax.XTick=[1 1.2 1.4];
ax.XTickLabel ={'All','M1','M2'};
ax.XTickLabelRotation=30;
xlim([0.8 1.8])
ylabel('Performance %')
set(gca,'fontsize',14,'fontweight','bold')

%% State for Performance
disp('%%% State for Performance %%%')

val_h=[]; val_h=performanceT;
Perfomance_all= round([nanmean(val_h),nanstd(val_h)/sqrt(length(val_h)),length(val_h)]);
disp (strcat('All Performance; ',' P',' =',num2str(Perfomance_all(1)),' + ',...
    num2str(Perfomance_all(2)),'%; n = ',num2str(Perfomance_all(3))));

val_h=[]; val_h=performanceT(ind_ses_m1);
Perfomance_m1= round([nanmean(val_h),nanstd(val_h)/sqrt(length(val_h)),length(val_h)]);
disp (strcat('Monkey 1 Performance; ',' P',' =',num2str(Perfomance_m1(1)),' + ',...
    num2str(Perfomance_m1(2)),'%; n = ',num2str(Perfomance_m1(3))));

val_h=[]; val_h=performanceT(ind_ses_m2);
Perfomance_m2= round([nanmean(val_h),nanstd(val_h)/sqrt(length(val_h)),length(val_h)]);
disp (strcat('Monkey 2 Performance; ',' P',' =',num2str(Perfomance_m2(1)),' + ',...
    num2str(Perfomance_m2(2)),'%; n = ',num2str(Perfomance_m2(3))));

%% Fig. 1b, Mean performance of individual monkeys
ind_ses=1;                                                               % All sessions  92
ind_ses_m1=setdiff(ind_ses,[52:92]);                                        % Sessions from 1-51 is for Monkey 1
ind_ses_m2=setdiff(ind_ses,[1:51]);                                         % Sessions from 52-92 is for Monkey 2

figure('name','Fig. 1b, Rezayat.E, et al')
ax=subplot(1,1,1);
hold on
bar([1],nanmean(performance_In(ind_ses)),0.2,'FaceColor',[0.1 0.1 0.1])
bar([1]+0.2,nanmean(performance_Out(ind_ses)),0.2,'FaceColor',[0.5 0.5 0.5])
% bar([1]+0.4,nanmean(performanceT(ind_ses_m2)),0.2,'FaceColor',[0.3 0.3 0.3])

val_h=[]; val_h=performance_In(ind_ses);
errorbar([1],nanmean(val_h),nanstd(val_h)/sqrt(length(val_h)),'color',[0.1 0.1 0.1],'linewidth',2)

val_h=[]; val_h=performance_Out(ind_ses);
errorbar([1]+0.2,nanmean(val_h),nanstd(val_h)/sqrt(length(val_h)),'color',[0.5 0.5 0.5],'linewidth',2)
% 
% val_h=[]; val_h=performanceT(ind_ses_m2);
% errorbar([1]+0.4,nanmean(val_h),nanstd(val_h)/sqrt(length(val_h)),'color',[0.3 0.3 0.3],'linewidth',2)

ylim([50 80])
ax.XTick=[1 1.2];
ax.XTickLabel ={'In','Out'};
ax.XTickLabelRotation=30;
xlim([0.8 1.5])
ylabel('Performance %')
set(gca,'fontsize',14,'fontweight','bold')


figure('name','Fig. 1b, Rezayat.E, et al')
ax=subplot(1,1,1);
hold on
bar([1],nanmean(performance_In(ind_ses)),0.2,'FaceColor',[0.1 0.1 0.1])
bar([1]+0.2,nanmean(performance_Out(ind_ses)),0.2,'FaceColor',[0.5 0.5 0.5])
% bar([1]+0.4,nanmean(performanceT(ind_ses_m2)),0.2,'FaceColor',[0.3 0.3 0.3])

val_h=[]; val_h=performance_In(ind_ses);
errorbar([1],nanmean(val_h),nanstd(val_h)/sqrt(length(val_h)),'color',[0.1 0.1 0.1],'linewidth',2)

val_h=[]; val_h=performance_Out(ind_ses);
errorbar([1]+0.2,nanmean(val_h),nanstd(val_h)/sqrt(length(val_h)),'color',[0.5 0.5 0.5],'linewidth',2)
% 
% val_h=[]; val_h=performanceT(ind_ses_m2);
% errorbar([1]+0.4,nanmean(val_h),nanstd(val_h)/sqrt(length(val_h)),'color',[0.3 0.3 0.3],'linewidth',2)

ylim([50 80])
ax.XTick=[1 1.2];
ax.XTickLabel ={'In','Out'};
ax.XTickLabelRotation=30;
xlim([0.8 1.5])
ylabel('Performance %')
set(gca,'fontsize',14,'fontweight','bold')



figure('name','Fig. S1a, Rezayat.E, et al')
p=signrank(performance_Out,performance_In);
n=length(performance_In);
hold on
scatter(performance_Out(ind_ses_m1),performance_In(ind_ses_m1),'b*')
scatter(performance_Out(ind_ses_m2),performance_In(ind_ses_m2),'gs')
plot([-100:250],[-100:250],'k')
ylabel('Preformance (%) [In]')
xlabel('Preformance (%) [Out]')
set(gca,'fontsize',14,'fontweight','bold')
 axis([48 95 48 95])

%% State for Performance
disp('%%% State for Performance %%%')

val_h=[]; val_h=performance_In;
Perfomance_all= round([nanmean(val_h),nanstd(val_h)/sqrt(length(val_h)),length(val_h)]);
disp (strcat('In Performance; ',' P',' =',num2str(Perfomance_all(1)),' + ',...
    num2str(Perfomance_all(2)),'%; n = ',num2str(Perfomance_all(3))));

val_h=[]; val_h=performance_Out(ind_ses);
Perfomance_m1= round([nanmean(val_h),nanstd(val_h)/sqrt(length(val_h)),length(val_h)]);
disp (strcat('Out Performance; ',' P',' =',num2str(Perfomance_m1(1)),' + ',...
    num2str(Perfomance_m1(2)),'%; n = ',num2str(Perfomance_m1(3))));

%% Fig. S1a, Reaction time for correct vs. wrong trials across sessions
ind_ses=1;                                                               % All sessions  86
ind_ses_m1=setdiff(ind_ses,[52:86]);                                        % Sessions from 1-51 is for Monkey 1
ind_ses_m2=setdiff(ind_ses,[1:51]);                                         % Sessions from 52-86 is for Monkey 2
p=signrank(RT_wr,RT_cr);
n=length(RT_cr);

figure('name','Fig. S1a, Rezayat.E, et al')

hold on
scatter(RT_wr(ind_ses_m1)/1000,RT_cr(ind_ses_m1)/1000,'b*')
scatter(RT_wr(ind_ses_m2)/1000,RT_cr(ind_ses_m2)/1000,'gs')
plot([-100:250],[-100:250],'k')
ylabel('RT (sec.) [Cr]')
xlabel('RT (sec.) [Wr]')
set(gca,'fontsize',14,'fontweight','bold')
axis([0.1 0.25 0.1 0.25])
%% State for reaction time
try
disp('%%% State for reaction time %%%')

x_h=RT_wr;y_h=RT_cr;
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('Reaction time, Cr vs. Wr, All; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));

x_h=RT_wr(ind_ses_m1);y_h=RT_cr(ind_ses_m1);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('Reaction time, Cr vs. Wr, M1; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));

x_h=RT_wr(ind_ses_m2);y_h=RT_cr(ind_ses_m2);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('Reaction time, Cr vs. Wr, M2; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));
catch
end
%% Fig. S1b, Variability in saccade landing points for correct vs. wrong trials across sessions
ind_ses=1;                                                               % All sessions  86
ind_ses_m1=setdiff(ind_ses,[52:86]);                                        % Sessions from 1-51 is for Monkey 1
ind_ses_m2=setdiff(ind_ses,[1:51]);                                         % Sessions from 52-86 is for Monkey 2

figure('name','Fig. S1b, Rezayat.E, et al')

hold on
scatter(SL_wr(ind_ses_m1),SL_cr(ind_ses_m1),'b*')
scatter(SL_wr(ind_ses_m2),SL_cr(ind_ses_m2),'gs')
plot([0:0.1:2.5],[0:0.1:2.5],'k')
axis([0 2.5 0 2.5])
ylabel('Sac.land.variation (dva) [Cr]')
xlabel('Sac.land.variation (dva) [Wr]')
set(gca,'fontsize',14,'fontweight','bold')

%% State for Variability in saccade landing points for correct vs. wrong trials across sessions
try
disp('%%% State for Variability in saccade landing points for correct vs. wrong trials across sessions %%%')
x_h=SL_wr;y_h=SL_cr;
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('Saccade Landing Points, Cr vs. Wr, All; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));

x_h=SL_wr(ind_ses_m1);y_h=SL_cr(ind_ses_m1);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('Saccade Landing Points, Cr vs. Wr, M1; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));

x_h=SL_wr(ind_ses_m2);y_h=SL_cr(ind_ses_m2);
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('Saccade Landing Points , Cr vs. Wr, M2; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));
catch
end
%% Fig. S1c Performance as a function of sample eccentricity; averaged across all trials for each session
ind_ses=1;                                                               % All sessions  86
ind_ses_m1=setdiff(ind_ses,[52:86]);                                        % Sessions from 1-51 is for Monkey 1
ind_ses_m2=setdiff(ind_ses,[1:51]);                                         % Sessions from 52-86 is for Monkey 2

sample_radius=sqrt(sum(sam_loc.^2,2));
figure('name','Fig. S1c, Rezayat.E, et al')
hold on

scatter(sample_radius(ind_ses_m1),performanceT(ind_ses_m1),'b*')
scatter(sample_radius(ind_ses_m2),performanceT(ind_ses_m2),'gs')
xlabel('Sample location (dva)')
ylabel('Performance (%)')
set(gca,'fontsize',14,'fontweight','bold')
xlim([0 13])
 set(gca,'fontsize',14,'fontweight','bold')

% [r,p]=correlation(x_h,y_h,1001);
% [r, m ,b]=regression(x_h,y_h,'one');
% n=length(x_h);
% xx=1:13;
% plot(xx,m*xx+b,'r')
% text( 1,60,strcat('r = ',num2str(r,2)),'fontsize',14,'fontweight','bold')
% text( 1,57,strcat('p = ',num2str(p,2)),'fontsize',14,'fontweight','bold')
% text( 1,54,strcat('n = ',num2str(length(x_h))),'fontsize',14,'fontweight','bold')

%% Fig. S1d Eye position during the task
try
ind_ses=1;%session_selection(1);                                               % All sessions  with enough number of wrong trials
ind_ses=setdiff(ind_ses,[87:92]);                                           % Eye data missed in 87:92
ind_ses_m1=setdiff(ind_ses,[52:86]);                                        % Sessions from 1-51 is for Monkey 1
ind_ses_m2=setdiff(ind_ses,[1:51]);                                         % Sessions from 52-86 is for Monkey 2

t_h=([1:2100]-500)/1000;
figure('name','Fig. S1d, Rezayat.E, et al')
ax=subplot(1,1,1);
hold on
eye_X_cr=squeeze((eye_pos_cr(:,:,1,:)));eye_X_wr=squeeze((eye_pos_wr(:,:,1,:)));
eye_Y_cr=squeeze((eye_pos_cr(:,:,2,:)));eye_Y_wr=squeeze((eye_pos_wr(:,:,2,:)));

niceplot(squeeze(eye_X_cr(ind_ses,4,1:end))+1,t_h,1,Color_loc(2,1),Color_loc(2,2),Color_loc(2,3))
niceplot(squeeze(eye_X_cr(ind_ses,1,1:end))+1,t_h,1,Color_loc(1,1),Color_loc(1,2),Color_loc(1,3))

niceplot(squeeze(eye_Y_cr(ind_ses,4,1:end))-1,t_h,1,Color_loc(2,1),Color_loc(2,2),Color_loc(2,3))
niceplot(squeeze(eye_Y_cr(ind_ses,1,1:end))-1,t_h,1,Color_loc(1,1),Color_loc(1,2),Color_loc(1,3))

niceplot(squeeze(eye_X_wr(ind_ses,1,1:end))-3,t_h,1,Color_bhv(2,1),Color_bhv(2,2),Color_bhv(2,3))
niceplot(squeeze(eye_X_cr(ind_ses,1,1:end))-3,t_h,1,Color_bhv(1,1),Color_bhv(1,2),Color_bhv(1,3))

niceplot(squeeze(eye_Y_wr(ind_ses,1,1:end))-5,t_h,1,Color_bhv(2,1),Color_bhv(2,2),Color_bhv(2,3))
niceplot(squeeze(eye_Y_cr(ind_ses,1,1:end))-5,t_h,1,Color_bhv(1,1),Color_bhv(1,2),Color_bhv(1,3))

xlabel('Time from sample onset (sec)')
ylabel('Eye poition (dva)')
xlim([-0.250 1.600])
line([0 0],ylim)
line([300 300]/1000,ylim)
line([1300 1300]/1000,ylim)
ylim([-6 2])
ax.XTick=[0 0.3 1.3];
ax.YTick=[-6:2];
ax.YTickLabel ={'','V','','H','','V','','H',''};
set(gca,'fontsize',14,'fontweight','bold')
text(0.7,1.7,'In','color',Color_loc(1,:),'fontsize',20)
text(0.7,0.5,'Out','color',Color_loc(2,:),'fontsize',20)
text(0.7,-2.2,'Cr','color',Color_bhv(1,:),'fontsize',20)
text(0.7,-3.5,'Wr','color',Color_bhv(2,:),'fontsize',20)
catch
end
%% Stat for Eye positions

try
    disp('%%% State for Eye positions for In vs. Out across sessions %%%')
x_h=squeeze(nanmean(eye_X_cr(ind_ses,4,t1_delay:t2_delay),3));
y_h=squeeze(nanmean(eye_X_cr(ind_ses,1,t1_delay:t2_delay),3));
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('Eye Position X, In vs. Out, All; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));

x_h=squeeze(nanmean(eye_Y_cr(ind_ses,4,t1_delay:t2_delay),3));
y_h=squeeze(nanmean(eye_Y_cr(ind_ses,1,t1_delay:t2_delay),3));
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('Eye Position Y, In vs. Out, All; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));

x_h=squeeze(nanmean(eye_X_wr(ind_ses,1,t1_delay:t2_delay),3));
y_h=squeeze(nanmean(eye_X_cr(ind_ses,1,t1_delay:t2_delay),3));
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('Eye Position X, Cr vs. Wr, ; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));

x_h=squeeze(nanmean(eye_Y_wr(ind_ses,1,t1_delay:t2_delay),3));
y_h=squeeze(nanmean(eye_Y_cr(ind_ses,1,t1_delay:t2_delay),3));
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('Eye Position Y, Cr vs. Wr, ; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));
catch
end
%% Fig. S1e Microsaccade rate over the timecourse of a trial for correct and wrong trials
try
ind_ses=1;%session_selection(1);                                               % All sessions  with enough number of wrong trials
ind_ses=setdiff(ind_ses,[87:92]);                                           % Eye data missed in 87:92
ind_ses_m1=setdiff(ind_ses,[52:86]);                                        % Sessions from 1-51 is for Monkey 1
ind_ses_m2=setdiff(ind_ses,[1:51]);                                         % Sessions from 52-86 is for Monkey 2
figure('name','Fig. S1d, Rezayat.E, et al')
ax=subplot(1,1,1);
hold on
niceplot(squeeze(micsac_wr(ind_ses,1,1:end)),t_h,10,Color_bhv(2,1),Color_bhv(2,2),Color_bhv(2,3))
niceplot(squeeze(micsac_cr(ind_ses,1,1:end)),t_h,10,Color_bhv(1,1),Color_bhv(1,2),Color_bhv(1,3))
xlabel('Time from sample onset (sec.)')
ylabel('Microsaccade rate(Hz)')
xlim([-0.25 1.6])
line([0 0],ylim)
line([0.3 0.3],ylim)
line([1.3 1.3],ylim)
ax.XTick=[0 0.3 1.3];
ax.YTick=[0 4];
set(gca,'fontsize',14,'fontweight','bold')
text(0.7,6,'Cr','color',Color_bhv(1,:),'fontsize',20)
text(0.7,2,'Wr','color',Color_bhv(2,:),'fontsize',20)
ylim([0 4])
catch
end
%% Stat for microsaccade
try
disp('%%% State for Microsaccade for Cr vs. Wr across sessions %%%')
x_h=squeeze(nanmean(micsac_wr(ind_ses,1,t1_delay:t2_delay),3));
y_h=squeeze(nanmean(micsac_cr(ind_ses,1,t1_delay:t2_delay),3));
p=signrank(x_h,y_h);
n=length(x_h);
disp (strcat('Microsaccade, Cr vs. Wr, All; ',' Diff',' =',num2str(nanmean(y_h)-nanmean(x_h)),' + ',...
    num2str(nanstd(y_h-x_h)/sqrt(n)),'; n = ',num2str(n),' p = ',num2str(p)));

catch
end
