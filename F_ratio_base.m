function F= F_ratio_base(pref,nonp,inter)
% This code calculates explained variance proposed by 
%Olejnik, S., & Algina, J. (2003). Generalized eta and omega 
%squared statistics: measures of effect size for some common
%research designs. Psychological methods, 8(4), 434.
if (pref==0); F=nan;return;end
if (nonp==0); F=nan;return;end
arr=[pref;nonp;inter];
labels_=[ones(length(pref),1); 2*ones(length(nonp),1); 3*ones(length(inter),1)];
[~,tb1]=anova1(arr,labels_,'off');


 F=(tb1{2,2}-(tb1{2,3}*tb1{3,4}))/(tb1{4,2}+tb1{3,4});
% F=tb1{2,5};

end