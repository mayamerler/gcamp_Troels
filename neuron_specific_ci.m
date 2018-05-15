
load('gcamp.mat');
load('rr5summary.mat');
% zone size
zs = 15;

%% before and after time zones

beforereward = mayabeforezones(rr5summary(1).m7677.rewarded_entries,zs);
afterreward = mayaafterzones(rr5summary(1).m7677.rewarded_entries,zs);
beforeunreward = mayabeforezones(rr5summary(1).m7677.unrewarded_entries,zs);
afterunreward = mayaafterzones(rr5summary(1).m7677.unrewarded_entries,zs);

%% auROC for just rewarded head entries

r_before=[];
for i = 1:size(beforereward,1)
    r_before(:,(i*zs-(zs-1)):(i*zs)) = gcamp(1).m7677.neuron(:,beforereward(i,1):beforereward(i,2));
    r_index_before((i*zs-(zs-1)):(i*zs)) = beforereward(i,1):beforereward(i,2);
end
r_before = r_before'; % makes plotting easier

r_after=[];
for i = 1:size(afterreward,1)
    r_after(:,(i*zs-(zs-1)):(i*zs)) = gcamp(1).m7677.neuron(:,afterreward(i,1):afterreward(i,2));
    r_index_after((i*zs-(zs-1)):(i*zs)) = afterreward(i,1):afterreward(i,2);
end
r_after = r_after';

auROCr = [];
auROCr = mayaauroc(gcamp(1).m7677.neuron,r_index_before,r_index_after);
[auROCr_inc, indexr] = sort(auROCr);
[ecdfr,xr] = ecdf(auROCr);
    
%% auROC for just unrewarded head entries

u_before = [];
for i = 1:size(beforeunreward,1)
    u_before(:,(i*zs-(zs-1)):(i*zs)) = gcamp(1).m7677.neuron(:,beforeunreward(i,1):beforeunreward(i,2));
    u_index_before((i*zs-(zs-1)):(i*zs)) = beforeunreward(i,1):beforeunreward(i,2);
end
u_before = u_before'; % makes plotting easier

u_after = [];
for i = 1:size(afterunreward,1)
    u_after(:,(i*zs-(zs-1)):(i*zs)) = gcamp(1).m7677.neuron(:,afterunreward(i,1):afterunreward(i,2));
    u_index_after((i*zs-(zs-1)):(i*zs)) = afterunreward(i,1):afterunreward(i,2);
end
u_after = u_after';

auROCu = [];
auROCu = mayaauroc(gcamp(1).m7677.neuron,u_index_before,u_index_after);
[auROCu_inc, indexu] = sort(auROCu);
auROCu_smooth = smooth(auROCu_inc,5);
[ecdfu,xu] = ecdf(auROCu);



%% randomizing times of events 
% clearvars auROCdr auROChr auROCpr auROCrr auROCur;

for f = 1:200 % Iterations for the randomization
    
    rand_reward = sort(randi(1500, size(rr5summary(1).m7677.rewarded_entries,1),1));
    rand_unreward = sort(randi(1500, size(rr5summary(1).m7677.unrewarded_entries,1),1));
    
% finding time zones before and after each randomized event
    
    beforerewardrand = mayabeforezones(rand_reward,zs);
    afterrewardrand = mayaafterzones(rand_reward,zs);
    
    beforeunrewardrand = mayabeforezones(rand_unreward,zs);
    afterunrewardrand = mayaafterzones(rand_unreward,zs);    


% auROC analysis on randomized rewarded entry times
    rr_before = [];
    rr_after = [];
    for i = 1:size(beforerewardrand,1)
        rr_before(:,(i*zs-(zs-1)):(i*zs)) = gcamp(1).m7677.neuron(:,beforerewardrand(i,1):beforerewardrand(i,2));
        rr_index_before((i*zs-(zs-1)):(i*zs)) = beforerewardrand(i,1):beforerewardrand(i,2);
    end
    rr_before = rr_before'; % makes plotting easier

    for i = 1:size(afterrewardrand,1)
        rr_after(:,(i*zs-(zs-1)):(i*zs)) = gcamp(1).m7677.neuron(:,afterrewardrand(i,1):afterrewardrand(i,2));
        rr_index_after((i*zs-(zs-1)):(i*zs)) = afterrewardrand(i,1):afterrewardrand(i,2);
    end
    rr_after = rr_after';
    
    auROCrr(f,:) = mayaauroc(gcamp(1).m7677.neuron,rr_index_before,rr_index_after);
    auROCrr_inc(f,:) = sort(auROCrr(f,:));

    
% auROC analysis on randomized unrewarded entry times
    ur_before = [];
    ur_after = [];
    for i = 1:size(beforeunrewardrand,1)
        ur_before(:,(i*zs-(zs-1)):(i*zs)) = gcamp(1).m7677.neuron(:,beforeunrewardrand(i,1):beforeunrewardrand(i,2));
        ur_index_before((i*zs-(zs-1)):(i*zs)) = beforeunrewardrand(i,1):beforeunrewardrand(i,2);
    end
    ur_before = ur_before'; % makes plotting easier
    
    for i = 1:size(afterrewardrand,1)
        ur_after(:,(i*zs-(zs-1)):(i*zs)) = gcamp(1).m7677.neuron(:,afterunrewardrand(i,1):afterunrewardrand(i,2));
        ur_index_after((i*zs-(zs-1)):(i*zs)) = afterunrewardrand(i,1):afterunrewardrand(i,2);
    end
    ur_after = ur_after';
    
    auROCur(f,:) = mayaauroc(gcamp(1).m7677.neuron,ur_index_before,ur_index_after);
    auROCur_inc(f,:) = sort(auROCur(f,:));
end
%%
clearvars auROCrr_CI tmp mus sigmas zscores
for i = 1:size(auROCrr,2)
    tmp = sort(auROCrr(:,i));
    auROCrr_CI(1,i) = tmp(5);
    auROCrr_CI(2,i) = tmp(195);
    [mus(i),sigmas(i)]= normfit(auROCrr_CI(:,i));
    zscores(i) = (auROCr(i)-mus(i))/sigmas(i);
end

z_pos = mean(zscores(find(zscores>0)));
z_neg = mean(zscores(find(zscores<0)));

%% distance from neuron-specific confidence interval
countneg = 0;
countpos = 0;
for i = 1:size(auROCr_inc,2)
    if auROCr(i) < auROCrr_CI(1,i)
        countneg = countneg+1;
        ci_distr_neg(countneg) = auROCrr_CI(i)-auROCr(i);
    end
    if auROCr(i) > auROCrr_CI(2,i)
        countpos = countpos + 1;
        ci_distr_pos(countpos) = auROCr_inc(i) - auROCrr_CI(i);
    end
end
if countpos + countneg == 0
ci_distr = 0; end
ci_outsider = (countneg + countpos)/i;
net_distr = sum(ci_distr_pos)+sum(ci_distr_neg);


