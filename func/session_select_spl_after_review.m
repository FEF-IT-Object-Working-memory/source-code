function ind_ses=session_select_spl_after_review(status)
%% session selection for comparing corect and Wrong
lfp_niose=[];
% lfp_niose=[20 21 23:37 41 60 88];
 lfp_niose= [20 21 23 24 25 27 28 29 31 32 33 34 35 37 60];
 
 ind_nonslective =[7  8  52  65  69 72  85 92];
lfp_niose=[lfp_niose ind_nonslective];
%   lfp_niose= [23:37 60];
if isempty(status)
    status =4;   % status = 1 All , status = 2 min3 wrong; status = 4 performance &3 min wrong
end
[conditions,performance]=session_stats_after_review();
switch status
    case 1   %% All
        ind_ses = [1:92];
          ind_ses=setdiff(ind_ses,lfp_niose);

    case 2  %% Number of minimum wrong trials >3
        [is_h ~] = find(sum(conditions.wr(:,[1 4]),2)<3);
        is_h = unique(is_h);
        ind_ses = [1:92];
        ind_ses = setdiff(ind_ses,is_h);
         ind_ses=setdiff(ind_ses,lfp_niose);
        % ind_ses = setdiff(ind_ses,17);
    case 3  %% 
        ind_ses = [1:92];
    case 4  %%  minimum wrong
  %       [is_h ~] = find(((performance(:,1)<60)|(performance(:,4)<60))|(performance(:,3)<60));
                
        
       [is_h ~] = find(sum(conditions.wr(:,[1 2 3]),2)<3);
        is_h = unique(is_h);
        ind_ses = [1:92];
        ind_ses = setdiff(ind_ses,is_h);
        
    case 5   %% All
        ind_ses = [1:92];
          ind_ses=setdiff(ind_ses,lfp_niose);

    case 6  %% Number of minimum wrong trials >3
        [is_h ~] = find(conditions.wr(:,1)<3);
        is_h = unique(is_h);
        ind_ses = [1:92];
        ind_ses = setdiff(ind_ses,is_h);
         ind_ses=setdiff(ind_ses,lfp_niose);
        % ind_ses = setdiff(ind_ses,17);
    case 7  %% 
        ind_ses = [1:92];
    case 8  %%  minimum wrong
  %       [is_h ~] = find(((performance(:,1)<60)|(performance(:,4)<60))|(performance(:,3)<60));
                
        
        [is_h ~] = find((conditions.wr(:,[1]))<3);
        is_h = unique(is_h);
        ind_ses = [1:92];
        ind_ses = setdiff(ind_ses,is_h);
        
    case 9  %% performance
        [is_h ~] = find((performance(:,[1:3])<60));
        is_h = unique(is_h);
        ind_ses = [1:92];
        ind_ses = setdiff(ind_ses,is_h);
    case 10
        [is_h ~] = find(sum(conditions.wr(:,[1,4]),2)<4);
        is_h = unique(is_h);
        ind_ses = [1:92];
        ind_ses = setdiff(ind_ses,is_h);
         ind_ses=setdiff(ind_ses,lfp_niose);
    case 11  %% performance
        [is_h ~] = find((performance(:,[1,4])<60));
        is_h = unique(is_h);
        ind_ses = [1:92];
        ind_ses = setdiff(ind_ses,is_h);
            
    case 12  %% performance & minimum wrong
         [is_h ~] = find(((performance(:,1)<50)));
%                 [is_h ~] = find(((performance(:,1)<60)));

        is_h = unique(is_h);
        [is_h2 ~] = find(conditions.wr(:,1)<3);
        is_h2 = unique(is_h2);
        ind_ses = [1:92];
        ind_ses = setdiff(ind_ses,is_h);
        ind_ses = setdiff(ind_ses,is_h2);
         ind_ses=setdiff(ind_ses,lfp_niose);
end

 %%ind_ses=setdiff(ind_ses,[23:37]); %%  This session It has alpha coherecny and power may from reward
%  ind_ses=setdiff(ind_ses,[3 8 11 15 16 22:36]); %%  This session It has alpha coherecny may from reward

 
 %ind_ses = setdiff(ind_ses,[1 2 32 33 60]); %%  This session It has no modulation in evoked response

 % ind_ses=setdiff(ind_ses,[38]);
  %    ind_ses=setdiff(ind_ses,[10 38 72 32]);

   % ind_ses=[10 38 72 32]
%  ind_ses=setdiff(ind_ses,[10 38 45 72]); %% This session low power beta and gamma

% ind_ses=setdiff(ind_ses,[10 15 24 32 34 35 38 45 50 62 72]); %% This session low power beta and gamma
%   ind_ses=[10 15 24 32 34 35 38 45 50 62 72];
% ind_ses=setdiff(ind_ses,[9:11 15 16 23:26 32:35 38:51 ]); %% This session low power beta and gamma

%ind_ses=setdiff(ind_ses,[1:8 38:45 ]); %% This session meial of IT
%ind_ses=[1:8 38:45 ];%% This session lateral of IT

% ind_ses=setdiff(ind_ses,[1:6 32:51 ]); %% This session anterior of IT
% ind_ses=[1:6 32:51 ];%% This session posterior of IT

%  ind_ses=setdiff(ind_ses,[7,8]);  % IT site is out  IT sit is M 15 is too close to frontal
%  ind_ses=setdiff(ind_ses,[1:22]);
%ind_ses=setdiff(ind_ses,[32:51 ]); %% This session lateral of FEF
%ind_ses=[32:51 ]; %% This sessions are medials of FEF
%[9:11 15 16 23:26 32:35 38:51 ];

