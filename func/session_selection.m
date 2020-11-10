function ind_ses=session_selection(status)
%  
%Data are reported from
%69 IT sites (35 M1, 34 M2) 
%92 FEF sites (51 M1, 41 M2). 
%LFP data of 15 IT sites (14 M1, 1 M2) was discarded prior to any
%analysis due to noise, 
%8 IT sites were excluded due to lack of object selective multiunit activity.
%When comparing correct and wrong trials, 
%certain sites and units were excluded due to low numbers of wrong trials (< 3 wrong trials)
% 1 [ALL], 2 [All Number of minimum wrong trials >3] 
it_lfp_niose=[]; it_lfp_niose= [20 21 23 24 25 27 28 29 31 32 33 34 35 37 60];
it_lfp_nonselective =[];  it_lfp_nonselective=[7  8  52  65  69 72  85 92  ];
%  it_lfp_nonselective=[ 7   8   52   65    69    72    86    87    88    89    91];
if isempty(status)
    status =1;
end
[conditions,performance]=session_stats();
switch status
    case 1                                                                  %% All
        ind_ses = [1:92];
        
    case 2                                                                  %% All Number of minimum wrong trials >3
        [is_h ~] = find(conditions.wr(:,1)<3);
        is_h = unique(is_h);
        ind_ses = [1:92];
        ind_ses = setdiff(ind_ses,is_h);
        
    case 3                                                                  %% All IT  LFP included
        ind_ses = [1:92];
        ind_ses=setdiff(ind_ses,it_lfp_niose);
        ind_ses=setdiff(ind_ses,it_lfp_nonselective);
        
    case 4                                                                  %% IT  LFP included, Number of minimum wrong trials >3
        [is_h ~] = find(conditions.wr(:,1)<3);
        is_h = unique(is_h);
        ind_ses = [1:92];
        ind_ses = setdiff(ind_ses,is_h);
        ind_ses=setdiff(ind_ses,it_lfp_niose);
        ind_ses=setdiff(ind_ses,it_lfp_nonselective);
% %         ind_ses = setdiff(ind_ses,59);
%   ind_ses = [ind_ses,59];
%     ind_ses = setdiff(ind_ses,[ 47]);
    case 5                                                                  %% Within FEF, Number of minimum wrong trials >3
        [is_h ~] = find(sum(conditions.wr(:,1:3),2)<3);
        is_h = unique(is_h);
        ind_ses = [1:92];
        ind_ses = setdiff(ind_ses,is_h);
        
    case 6                                                                  %% Within IT no LFP, Number of minimum wrong trials >3
        [is_h ~] = find(sum(conditions.wr(:,[1,4]),2)<3);
        is_h = unique(is_h);
        ind_ses = [1:92];
        ind_ses = setdiff(ind_ses,is_h);
        
    case 7                                                                  %% Within IT LFP, Number of minimum wrong trials >3
        [is_h ~] = find(sum(conditions.wr(:,[1,4]),2)<3);
        is_h = unique(is_h);
        ind_ses = [1:92];
        ind_ses = setdiff(ind_ses,is_h);
        ind_ses=setdiff(ind_ses,it_lfp_niose);
        ind_ses=setdiff(ind_ses,it_lfp_nonselective);
        
end

