function [conditions,performance]=session_stats_after_review()
  load('F:\Data\Reaserch\Thesis\Data Analysis\Data\gen_session_inf_500ms_baseline_full.mat')
% load('F:\Data\Reaserch\Thesis\Data Analysis Reconfiguration\Data\gen_session_inf.mat')

conditions = [];
for ss =1:92
    obj_ix_h = [gen_session_inf(ss).preff_obj,gen_session_inf(ss).npreff_obj];
    obj_ix_h = [obj_ix_h(1) setdiff([1:3],obj_ix_h) obj_ix_h(2)];
    for condi =1:6
        if condi<4; obj_h = obj_ix_h(condi);else obj_h=obj_ix_h(condi-3)+3;end
        ix_h = [];
%        if ss<74
        
        ix_h = ismember(gen_session_inf(ss).condition{:},obj_h)&(gen_session_inf(ss).bhv{:}==1)&(~gen_session_inf(ss).FEF_artifact)&(~gen_session_inf(ss).IT_artifact);
%        else
%                  ix_h = ismember(gen_session_inf(ss).condition{:},obj_h)&(gen_session_inf(ss).bhv{:}==1)&(~gen_session_inf(ss).FEF_artifact);
%        end  
        conditions.cr(ss,condi) = sum(ix_h);
        ix_h = [];        
        
        ix_h = ismember(gen_session_inf(ss).condition{:},obj_h)&(gen_session_inf(ss).bhv{:}==0)&(~gen_session_inf(ss).FEF_artifact)&(~gen_session_inf(ss).IT_artifact);
         
        conditions.wr(ss,condi) = sum(ix_h);
    end
end
performance = (conditions.cr)./(conditions.cr+conditions.wr)*100;

end