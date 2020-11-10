function out = gen_fx_get_svm(grp0,I0,rate,rep)
Pt = [];
Tu = [];
thr = 10000;
for cnt = 1:rep
    if length(grp0) > thr
        ix = randperm(length(grp0),thr);
        grp = grp0(ix); I = I0(ix,:);
    else
        grp = grp0; I = I0;
    end
    
    [test train] = gen_fx_get_equal_part(grp,rate);
    if sum(test)>5 && sum(train) > 5
        cls = gen_fx_MC_SVM(I(test,:),I(train,:),grp(train));
        Pt = [Pt; sum(cls == grp(test))/sum(test)];
        Ct(:,:,cnt)= confusionmat(cls,grp(test));
        Tu = [Tu diag(Ct(:,:,cnt)) ./ (sum(Ct(:,:,cnt)))'];
    else
        Pt = [Pt; nan];
        Ct(:,:,cnt)= nan;
        Tu = [Tu nan];
    end
end
out.C = Ct; out.pt = Pt; out.tu = Tu;
end