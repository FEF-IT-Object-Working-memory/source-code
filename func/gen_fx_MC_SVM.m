function class = gen_fx_MC_SVM(sample,train,grp)

classNo = double(unique(grp));
cls_nu = length(classNo);
%r = combntns(1:cls_nu,2);
r = nchoosek(1:cls_nu,2);
vote = [];
for i =1 : length(r(:,1))
    ix = (grp == classNo(r(i,1))) | (grp == classNo(r(i,2)));
    svmStruct = fitcsvm(train(ix,:),grp(ix));
    vote = [vote predict(svmStruct,sample)];
end
class = mode(vote,2);

% % Mdl = fitcecoc(train,grp);
% % class = predict(Mdl,sample);
