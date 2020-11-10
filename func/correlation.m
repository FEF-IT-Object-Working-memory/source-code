function [r,p]=correlation(x,y,P)

if size(x)&size(y)
    val= (~isnan(x))&(~isnan(y));
    x=x(val);y=y(val);
end
if size(x)&size(y)
   % r_ = []; p_ = [];
%    [r_,p_]= corrcoef(x,y,'rows','pairwise');
%     corr_=r_(1,2);p=p_(1,2);
    [r_,p_]= corr(x,y,'type','kendal');
    corr_=r_;p=p_;
    
    x_=x;corr_p=0;
    for o=1:P
        ind_ = [];ind_=randperm(size(x,1));
        x_=x(ind_);
         ind_ = [];ind_=randperm(size(y,1));
        y_=y(ind_);
%         r_ = []; r_=corrcoef(x_,y_,'rows','pairwise');
%         corr_p(o)=r_(1,2);
           r_ = []; r_=corr(x_,y_,'type','kendal');
        corr_p(o)=r_;
    end
    r=corr_-mean(corr_p);
    if corr_>0
    p_val=size(find(corr_p>=corr_))/size(corr_p);
    else
    p_val=size(find(corr_p<=corr_))/size(corr_p);  
    end
    p=p_val;
else
    r= nan; p = nan;
end

