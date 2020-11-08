function resp=Preprocessing(resp)

%% Reformation of conditions [Combine conditions]
gen_cond=resp.condition;
for i=1:size( gen_cond,1);
    if(gen_cond(i)<5)
        gen_cond_(i)=1;
    elseif (gen_cond(i)<9)
        gen_cond_(i)=2;
    elseif (gen_cond(i)<13)
        gen_cond_(i)=3;
    elseif (gen_cond(i)<17)
        gen_cond_(i)=4;
    elseif (gen_cond(i)<21)
        gen_cond_(i)=5;
    elseif (gen_cond(i)<25)
        gen_cond_(i)=6;
    end
end
resp.condition=gen_cond_;

%% Preferred Object

it_units=resp(1).IT;
mu=it_units{1};
for uini=1:size(it_units,2)
    obj_ind=find_preferance(it_units{uini},resp(1).condition);
    pref_h(uini)=obj_ind.pref;
    npref_h(uini)=obj_ind.npref;
    if uini>1
        mu=mu|it_units{uini};
    end
end
resp(1).pref_obj_su=pref_h;
resp(1).npref_obj_su=npref_h;

obj_ind=find_preferance(mu,resp(1).condition);

resp(1).pref_obj_mu=obj_ind.pref;
resp(1).npref_obj_mu=obj_ind.npref;

%% LFP artifact
resp=LFP_Artifact_remove(resp);

end