%% Loads data for mouse 7677 over 5 days

files{1} = 'D:\7677 gcamp analysis\7677 analysis\Workspace 7677 RR5 Day1 calculated.mat';
files{2} = 'D:\7677 gcamp analysis\7677 analysis\Workspace 7677 RR5 Day2 calculated.mat';
files{3} = 'D:\7677 gcamp analysis\7677 analysis\Workspace 7677 RR5 Day3 calculated.mat';
files{4} = 'D:\7677 gcamp analysis\7677 analysis\Workspace 7677 RR5 Day4 calculated.mat';
files{5} = 'D:\7677 gcamp analysis\7677 analysis\Workspace 7677 RR5 Day5.mat';

% List of figures generated:

% Let's fill this in...
%
%

% Other needed edits:
% (1) Get rid of cell arrays in z structure (just make them regular arrays)
% (2) Simplify code by substituting mayaauroc function where appropriate
% (3)  
%
% Use allfigures = 0 to suppress figures (only analyze data)

% List of dependencies:
%
% Other functions needed for this to run, e.g., mayabeforezones (which
% themselves need to be a bit more commented)

% Parameters for analysis
day=5; % Set day to a value between 1 and 5, run the whole script to generate figures
zs = 15; % sets size of time zone in samples (sampling rate for this data is 5 Hz)

load(files{day});
load('E:\rr5data.mat'); % This is to get the RR5 data without Excel 

widegraph = [57.8000  257.8000  990.4000  504.0000];

%% Initializing 

set(0,'defaultlinelinewidth',1.5);
set(0,'defaultfigurecolor','w');

allfigures=1;
clearvars auROC* z

%% average firing for C_raw

% Unsmoothed, raw average
avgC = mean(neuron.C_raw,1);
if allfigures; figure(1); plot(avgC); end

% A smoothed version
smooth_avgC = smooth(avgC,5);
if allfigures; figure(2); plot(smooth_avgC); end

%% compare binned data to smooth data

bindata = @(data,binsize) sum(reshape(data,binsize,numel(data)/binsize),1);
binned_avgC = bindata(avgC,10);
if allfigures; figure(3); plot(binned_avgC); end;

%% imports time data from excel file (rows and columns may vary)

rr5_7677 = rr5data(day).data([3:end], [8:9:26]);
rr5_7677(find(rr5_7677(:,1)==0),:)=[]; % Eliminate random zeroes

rr5head = rr5_7677(:,1)*5; %converts seconds to frames and caps at 1500 frames
rr5dipper = rr5_7677(:,2)*5;
rr5press = rr5_7677(:,3)*5;
z.head_entries{3} = rr5head(find(rr5head<1500));
z.dip_present{3} = rr5dipper(find(rr5dipper<1500));
z.rlev_press{3} = rr5press(find(rr5press<1500));

%% separate rewarded from non-rewarded head entries

rcounter = 1;
ucounter = 1;
for i = 1:size(z.head_entries{3},1)
    if sum((z.head_entries{3}(i) > z.dip_present{3}) & (z.head_entries{3}(i) < (z.dip_present{3}+25)))>0
        z.rewarded_entries{3}(rcounter,1) = z.head_entries{3}(i);
        rcounter = rcounter + 1;
    else
        z.unrewarded_entries{3}(ucounter,1) = z.head_entries{3}(i);
        ucounter = ucounter + 1;
    end
end

%% finding time zones before and after each event

for i = 3
    
    z.beforepress{i} = mayabeforezones(z.rlev_press{i},zs);
    z.afterpress{i} = mayaafterzones(z.rlev_press{i},zs);
    z.beforedip{i} = mayabeforezones(z.dip_present{i},zs);
    z.afterdip{i} = mayaafterzones(z.dip_present{i},zs);
    z.beforehead{i} = mayabeforezones(z.head_entries{i},zs);
    z.afterhead{i} = mayaafterzones(z.head_entries{i},zs);
    z.beforereward{i} = mayabeforezones(z.rewarded_entries{i},zs);
    z.afterreward{i} = mayaafterzones(z.rewarded_entries{i},zs);
    z.beforeunreward{i} = mayabeforezones(z.unrewarded_entries{i},zs);
    z.afterunreward{i} = mayaafterzones(z.unrewarded_entries{i},zs);
    
end

%% auROC for lever presses
% grouping all data from 'before' zones and all data from 'after' zones
p_before = [];
for i = 1:size(z.beforepress{3},1)
    p_before(:,(i*zs-(zs-1)):(i*zs)) = neuron.C_raw(:,z.beforepress{3}(i,1):z.beforepress{3}(i,2));
end
p_before = p_before'; % makes plotting easier

p_after = [];
for i = 1:size(z.afterpress{3},1)
        p_after(:,(i*zs-(zs-1)):(i*zs)) = neuron.C_raw(:,z.afterpress{3}(i,1):z.afterpress{3}(i,2));
end
p_after = p_after';

% find max, min, and increment that gives 100 thresholds between the two
highestpoint = max([p_before(:);p_after(:)]);
lowestpoint = min([p_before(:);p_after(:)]);
increment = (highestpoint - lowestpoint)/100; 
auROC = [];

% loop through neurons
for j = 1:size(neuron.C_raw,1)
    
    combinedata = [p_before(:,j);p_after(:,j)];
    codes = [zeros(numel(p_before(:,j)),1);ones(numel(p_after(:,j)),1)]; 
    
    %loop through thresholds
    for i = 1:100
        threshold = i*increment + lowestpoint;
        result = gt(combinedata,threshold);
        %calculate true positives and false positives
        tp(i)=sum((result==1)&(codes==1))/numel(p_after(:,j));
        fp(i)=sum((result==1)&(codes==0))/numel(p_before(:,j));
    end
    
    % Calculate the Riemann sum underneath the curve (AUC)
    [a,b] = sort(fp); 
    auROC(j) = trapz(sort(fp),tp(b));
    
end

% smooth and plot auROC
[auROC_inc, index] = sort(auROC);
auROC_smooth = smooth(auROC_inc,5);
if allfigures; figure; plot(auROC_inc, '.-');  title('auROC lever presses'); end

%% repeat for dipper presentations
d_before = [];
for i = 1:size(z.beforedip{3},1)
    d_before(:,(i*zs-(zs-1)):(i*zs)) = neuron.C_raw(:,z.beforedip{3}(i,1):z.beforedip{3}(i,2));
end
d_before = d_before'; % makes plotting easier
d_after = [];
for i = 1:size(z.afterdip{3},1)
        d_after(:,(i*zs-(zs-1)):(i*zs)) = neuron.C_raw(:,z.afterdip{3}(i,1):z.afterdip{3}(i,2));
end
d_after = d_after';

% find max, min, and increment that gives 100 thresholds between the two
highestpointd = max([d_before(:);d_after(:)]);
lowestpointd = min([d_before(:);d_after(:)]);
incrementd = (highestpointd - lowestpointd)/100; 
auROCd = [];

% loop through neurons
for j = 1:size(neuron.C_raw,1)
    
    combinedatad = [d_before(:,j);d_after(:,j)];
    codesd = [zeros(numel(d_before(:,j)),1);ones(numel(d_after(:,j)),1)]; 
    
    %loop through thresholds
    for i = 1:100
        threshold = i*incrementd + lowestpointd;
        result = gt(combinedatad,threshold);
        %calculate true positives and false positives
        tpd(i)=sum((result==1)&(codesd==1))/numel(d_after(:,j));
        fpd(i)=sum((result==1)&(codesd==0))/numel(d_before(:,j));
    end
    
    % Calculate the Riemann sum underneath the curve (AUC)
    [a,b] = sort(fpd); 
    auROCd(j) = trapz(sort(fpd),tpd(b));
    
end

% smooth and plot auROC
[auROCd_inc, indexd] = sort(auROCd);
auROCd_smooth = smooth(auROCd_inc,5);
if allfigures; figure; plot(auROCd_inc, '.-');  title('auROC dipper presentations'); end

%% repeat for head entries

h_before = [];
h_after = [];
for i = 1:size(z.beforehead{3},1)
    h_before(:,(i*zs-(zs-1)):(i*zs)) = neuron.C_raw(:,z.beforehead{3}(i,1):z.beforehead{3}(i,2));
end
h_before = h_before'; % makes plotting easier

for i = 1:size(z.afterhead{3},1)
        h_after(:,(i*zs-(zs-1)):(i*zs)) = neuron.C_raw(:,z.afterhead{3}(i,1):z.afterhead{3}(i,2));
end
h_after = h_after';

% find max, min, and increment that gives 100 thresholds between the two
highestpointh = max([h_before(:);h_after(:)]);
lowestpointh = min([h_before(:);h_after(:)]);
incrementh = (highestpointh - lowestpointh)/100; 
auROCh = [];

% loop through neurons
for j = 1:size(neuron.C_raw,1)
    
    combinedatah = [h_before(:,j);h_after(:,j)];
    codesh = [zeros(numel(h_before(:,j)),1);ones(numel(h_after(:,j)),1)]; 
    
    %loop through thresholds
    for i = 1:100
        threshold = i*incrementh + lowestpointh;
        result = gt(combinedatah,threshold);
        %calculate true positives and false positives
        tph(i)=sum((result==1)&(codesh==1))/numel(h_after(:,j));
        fph(i)=sum((result==1)&(codesh==0))/numel(h_before(:,j));
    end
    
    % Calculate the Riemann sum underneath the curve (AUC)
    [a,b] = sort(fph); 
    auROCh(j) = trapz(sort(fph),tph(b));
    
end
auROCh_inc = sort(auROCh);

%% auROC for just rewarded head entries

r_before=[];
for i = 1:size(z.beforereward{3},1)
    r_before(:,(i*zs-(zs-1)):(i*zs)) = neuron.C_raw(:,z.beforereward{3}(i,1):z.beforereward{3}(i,2));
end
r_before = r_before'; % makes plotting easier

r_after=[];
for i = 1:size(z.afterreward{3},1)
        r_after(:,(i*zs-(zs-1)):(i*zs)) = neuron.C_raw(:,z.afterreward{3}(i,1):z.afterreward{3}(i,2));
end
r_after = r_after';

% find max, min, and increment that gives 100 thresholds between the two
highestpointr = max([r_before(:);r_after(:)]);
lowestpointr = min([r_before(:);r_after(:)]);
incrementr = (highestpointr - lowestpointr)/100; 
auROCr = [];

% loop through neurons
for j = 1:size(neuron.C_raw,1)
    
    combinedatar = [r_before(:,j);r_after(:,j)];
    codesr = [zeros(numel(r_before(:,j)),1);ones(numel(r_after(:,j)),1)]; 
    
    %loop through thresholds
    for i = 1:100
        threshold = i*incrementr + lowestpointr;
        result = gt(combinedatar,threshold);
        %calculate true positives and false positives
        tpr(i)=sum((result==1)&(codesr==1))/numel(r_after(:,j));
        fpr(i)=sum((result==1)&(codesr==0))/numel(r_before(:,j));
    end
    
    % Calculate the Riemann sum underneath the curve (AUC)
    [a,b] = sort(fpr); 
    auROCr(j) = trapz(sort(fpr),tpr(b));
    
end

% smooth and plot auROC
[auROCr_inc, indexr] = sort(auROCr);
auROCr_smooth = smooth(auROCr_inc,5);

%% auROC for just unrewarded head entries

u_before = [];
for i = 1:size(z.beforeunreward{3},1)
    u_before(:,(i*zs-(zs-1)):(i*zs)) = neuron.C_raw(:,z.beforeunreward{3}(i,1):z.beforeunreward{3}(i,2));
end
u_before = u_before'; % makes plotting easier

u_after = [];
for i = 1:size(z.afterunreward{3},1)
        u_after(:,(i*zs-(zs-1)):(i*zs)) = neuron.C_raw(:,z.afterunreward{3}(i,1):z.afterunreward{3}(i,2));
end
u_after = u_after';

% find max, min, and increment that gives 100 thresholds between the two
highestpointu = max([u_before(:);u_after(:)]);
lowestpointu = min([u_before(:);u_after(:)]);
incrementu = (highestpointu - lowestpointu)/100; 
auROCu = [];

% loop through neurons
for j = 1:size(neuron.C_raw,1)
    
    combinedatau = [u_before(:,j);u_after(:,j)];
    codesu = [zeros(numel(u_before(:,j)),1);ones(numel(u_after(:,j)),1)]; 
    
    %loop through thresholds
    for i = 1:100
        threshold = i*incrementu + lowestpointu;
        result = gt(combinedatau,threshold);
        %calculate true positives and false positives
        tpu(i)=sum((result==1)&(codesu==1))/numel(u_after(:,j));
        fpu(i)=sum((result==1)&(codesu==0))/numel(u_before(:,j));
    end
    
    % Calculate the Riemann sum underneath the curve (AUC)
    [a,b] = sort(fpu); 
    auROCu(j) = trapz(sort(fpu),tpu(b));
    
end

% smooth and plot auROC
[auROCu_inc, indexu] = sort(auROCu);
auROCu_smooth = smooth(auROCu_inc,5);
figure; plot(auROCu_inc, '.-','markersize',12);  title('auROC unrewarded/rewarded head entries');
hold on; plot(auROCr_inc, 'r.-','markersize',12);
%% plot auROC with separate rewarded & unrewarded entries

ylims = [0,1];

figure; subplot(2,2,1); plot(auROC_inc, 'o-'); title('lever presses'); ylim(ylims);
subplot(2,2,2); plot(auROCd_inc, 'o-');  title('dipper presentations'); ylim(ylims);
subplot(2,2,3); plot(auROCh_inc, 'o-'); title('all head entries'); ylim(ylims); 
subplot(2,2,4); plot(auROCu_inc, 'ro-'); title('head entries by reward'); ylim(ylims);hold on;
subplot(2,2,4); plot(auROCr_inc, 'go-'); ylim(ylims); 
%set(gcf,'position',boxgraph)

saveas(gcf,'auROC_by_event.png');

%% averaged data with stimulus info

figure; subplot(3,1,1); plot(smooth_avgC, 'k'); title('averaged data with head entries'); hold on; 

for h = 1:size(z.head_entries{3},1)
    plot([z.head_entries{3}(h,1),z.head_entries{3}(h,1)],[-0.2,0.7], '--r','linewidth',0.8);
end

subplot(3,1,2); plot(smooth_avgC, 'k'); title('averaged data with lever presses'); hold on;
for p = 1:size(z.rlev_press{3},1)
    plot([z.rlev_press{3}(p,1),z.rlev_press{3}(p,1)],[-0.2,0.7],'--g','linewidth',0.8);
end

subplot(3,1,3); plot(smooth_avgC, 'k'); title('averaged data with dipper presentations'); hold on;
for i = 1:size(z.dip_present{3},1)
    plot([z.dip_present{3}(i,1),z.dip_present{3}(i,1)],[-0.2,0.7], '--b','linewidth',0.8);
end 

%set(gcf,'position',widegraph)
saveas(gcf,'event_trigged_average.png');
%% averaged data with head entries & dipper presentations
figure; plot(smooth_avgC, 'k'); title('averaged data with head entries (r) & dipper presentations (b)'); hold on;
for h = 1:size(z.head_entries{3},1)
    plot([z.head_entries{3}(h,1),z.head_entries{3}(h,1)],[-0.2,0.7], '--r','linewidth',0.8);
end
hold on;
for i = 1:size(z.dip_present{3},1)
    plot([z.dip_present{3}(i,1),z.dip_present{3}(i,1)],[-0.2,0.7], '--b','linewidth',0.8);
end 
%set(gcf,'position',widegraph)
saveas(gcf,'dipper_trigged_average.png');

%% maintaining same order of neurons over all plots (indexed by dipper presentation order)

figure; subplot(1,4,2); plot(auROC(indexd), 'o'); title('lever presses'); ylim([0.1,0.9]);
subplot(1,4,1); plot(auROCd(indexd), 'ro');  title('dipper presentations'); ylim([0.1,0.9]);
ylabel(sprintf('Day %1.0f auROC per neuron\n (in order of increasing value for reward)\n with before and after comparison window of %1.0f s',day,zs/5))
subplot(1,4,3); plot(auROCr(indexd), 'o'); title('rewarded head entries'); ylim([0.1,0.9]);
subplot(1,4,4); plot(auROCu(indexd), 'o'); title('unrewarded head entries'); ylim([0.1,0.9]);

set(gcf,'position',widegraph)
saveas(gcf,'auROC_indexed_to_dippervalues.png');

%% maintaining same order of neurons over all plots (indexed by rewarded entry order)

figure; subplot(1,4,2); plot(auROC(indexr), 'o'); title('lever presses'); ylim([0.1,0.9]);
subplot(1,4,3); plot(auROCd(indexr), 'o');  title('dipper presentations'); ylim([0.1,0.9]);
subplot(1,4,1); plot(auROCr(indexr), 'ro'); title('rewarded head entries'); ylim([0.1,0.9]);
ylabel(sprintf('Day %1.0f auROC per neuron\n (in order of increasing value for reward)\n with before and after comparison window of %1.0f s',day,zs/5))
subplot(1,4,4); plot(auROCu(indexr), 'o'); title('unrewarded head entries'); ylim([0.1,0.9]);

set(gcf,'position',widegraph)
saveas(gcf,'auROC_indexed_to_reward.png');

%% maintaining same order of neurons over all plots (indexed by lever press order)
figure; subplot(1,4,1); plot(auROC(index), 'ro'); title('lever presses'); ylim([0.1,0.9]);
ylabel(sprintf('Day %1.0f auROC per neuron\n (in order of increasing value for reward)\n with before and after comparison window of %1.0f s',day,zs/5))
subplot(1,4,2); plot(auROCd(index), 'o');  title('dipper presentations'); ylim([0.1,0.9]);
subplot(1,4,3); plot(auROCr(index), 'o'); title('rewarded head entries'); ylim([0.1,0.9]);
subplot(1,4,4); plot(auROCu(index), 'o'); title('unrewarded head entries'); ylim([0.1,0.9]);
set(gcf,'position',widegraph)
saveas(gcf,'auROC_indexed_to_leverpress.png');

%% maintaining same order of neurons over all plots (indexed by unrewarded entry order)
figure; subplot(1,4,2); plot(auROC(indexu), 'o'); title('lever presses'); ylim([0.1,0.9]);
subplot(1,4,4); plot(auROCd(indexu), 'o');  title('dipper presentations'); ylim([0.1,0.9]);
subplot(1,4,3); plot(auROCr(indexu), 'o'); title('rewarded head entries'); ylim([0.1,0.9]);
subplot(1,4,1); plot(auROCu(indexu), 'ro'); title('unrewarded head entries'); ylim([0.1,0.9]);
ylabel(sprintf('Day %1.0f auROC per neuron\n (in order of increasing value for reward)\n with before and after comparison window of %1.0f s',day,zs/5))
set(gcf,'position',widegraph)
saveas(gcf,'auROC_indexed_to_reward.png');

%% correlation coefficients between groups
figure; imagesc(corrcoef([auROC(indexd);auROCd(indexd);auROCr(indexd);auROCu(indexd)]'));
xticks([1:4]); yticks([1:4]); xticklabels({'lever','dipper','rewarded','unrewarded'});
yticklabels({'lever','dipper','rewarded','unrewarded'}); 
title(sprintf('Correlation coefficients (neuron auROC x neuron auROC)\nfor different periods of day %1.0f\n with before/after comparison window of %1.0f s',day,zs/5));
colorbar();

%% Activity of two neurons with auROC over 0.7 in dipper presentation data
for j = 1:size(neuron.C_raw,1)
    smoothed_rawdata(j,:) = smooth(neuron.C_raw(j,:),5);
end

cutoff = mean(auROCd)+std(auROCd);
figure; subplot(2,2,1); plot(smoothed_rawdata(auROCd>cutoff,:)'); hold on; 
for i = 1:size(z.dip_present{3},1)
    plot([z.dip_present{3}(i,1),z.dip_present{3}(i,1)],[-2,7], 'k');
end
title(sprintf('%1.0f neurons with auROC>%1.2f (locked to dipper)',sum(auROCd>cutoff),cutoff));

cutoff = mean(auROCd)-std(auROCd);
subplot(2,2,2); plot(smoothed_rawdata(auROCd<cutoff,:)'); hold on; 
for i = 1:size(z.dip_present{3},1)
    plot([z.dip_present{3}(i,1),z.dip_present{3}(i,1)],[-2,7], 'k');
end
title(sprintf('%1.0f neurons with auROC<%1.2f (locked to dipper)',sum(auROCd<cutoff),cutoff));

cutoff = mean(auROCd)+std(auROCd);
subplot(2,2,3); plot(mean(smoothed_rawdata(auROCd>cutoff,:)',2)); hold on; 
for i = 1:size(z.dip_present{3},1)
    plot([z.dip_present{3}(i,1),z.dip_present{3}(i,1)],[-2,7], 'k');
end
ylim([-1,3]);
title(sprintf('%1.0f neurons with auROC>%1.2f (locked to dipper)',sum(auROCd>cutoff),cutoff));

cutoff = mean(auROCd)-std(auROCd);
subplot(2,2,4); plot(mean(smoothed_rawdata(auROCd<cutoff,:)',2)); hold on; 
for i = 1:size(z.dip_present{3},1)
    plot([z.dip_present{3}(i,1),z.dip_present{3}(i,1)],[-2,7], 'k');
end
ylim([-1,3]);
title(sprintf('%1.0f neurons with auROC<%1.2f (locked to dipper)',sum(auROCd<cutoff),cutoff));
set(gcf,'position',widegraph)
saveas(gcf,'averages_by_auROC.png');

%% randomizing times of events 
% clearvars auROCdr auROChr auROCpr auROCrr auROCur;

for f = 1:200 % Iterations for the randomization
    
    rand_dip = sort(randi(1500, size(z.dip_present{3},1),1));
    rand_press = sort(randi(1500, size(z.rlev_press{3},1),1));
    rand_head = sort(randi(1500, size(z.head_entries{3},1),1));
    rand_reward = sort(randi(1500, size(z.rewarded_entries{3},1),1));
    rand_unreward = sort(randi(1500, size(z.unrewarded_entries{3},1),1));
    
% finding time zones before and after each randomized event
    
    z.beforepressrand{3} = mayabeforezones(rand_press,zs);
    z.afterpressrand{3} = mayaafterzones(rand_press,zs);
    
    z.beforediprand{3} = mayabeforezones(rand_dip,zs);
    z.afterdiprand{3} = mayaafterzones(rand_dip,zs);
    
    z.beforeheadrand{3} = mayabeforezones(rand_head,zs);
    z.afterheadrand{3} = mayaafterzones(rand_head,zs);
    
    z.beforerewardrand{3} = mayabeforezones(rand_reward,zs);
    z.afterrewardrand{3} = mayaafterzones(rand_reward,zs);
    
    z.beforeunrewardrand{3} = mayabeforezones(rand_unreward,zs);
    z.afterunrewardrand{3} = mayaafterzones(rand_unreward,zs);    

% repeat auROC analysis on randomized dipper presentation times
    dr_before = [];
    dr_after = [];
    for i = 1:size(z.beforediprand{3},1)
        dr_before(:,(i*zs-(zs-1)):(i*zs)) = neuron.C_raw(:,z.beforediprand{3}(i,1):z.beforediprand{3}(i,2));
    end
    dr_before = dr_before'; % makes plotting easier

    for i = 1:size(z.afterdiprand{3},1)
        dr_after(:,(i*zs-(zs-1)):(i*zs)) = neuron.C_raw(:,z.afterdiprand{3}(i,1):z.afterdiprand{3}(i,2));
    end
    dr_after = dr_after';

% find max, min, and increment that gives 100 thresholds between the two
    highestpointdr = max([dr_before(:);dr_after(:)]);
    lowestpointdr = min([dr_before(:);dr_after(:)]);
    incrementdr = (highestpointdr - lowestpointdr)/100; 
    % loop through neurons
    for j = 1:size(neuron.C_raw,1)
        combinedatadr = [dr_before(:,j);dr_after(:,j)];
        codesdr = [zeros(numel(dr_before(:,j)),1);ones(numel(dr_after(:,j)),1)]; 
    
    %loop through thresholds
        for i = 1:100
            threshold = i*incrementdr + lowestpointdr;
            result = gt(combinedatadr,threshold);
            %calculate true positives and false positives
            tpdr(i)=sum((result==1)&(codesdr==1))/numel(dr_after(:,j));
            fpdr(i)=sum((result==1)&(codesdr==0))/numel(dr_before(:,j));
        end
    
        % Calculate the Riemann sum underneath the curve (AUC)
        [a,b] = sort(fpdr); 
        auROCdr(f,j) = trapz(sort(fpdr),tpdr(b));
    
    end

    % smooth and plot auROC
    auROCdr_inc(f,:) = sort(auROCdr(f,:));
    %figure; plot(auROCdr_inc, '.-');  title('auROC with randomized dipper presentations');

    
% repeat auROC analysis on randomized lever press times
    pr_before = [];
    pr_after = [];
    for i = 1:size(z.beforepressrand{3},1)
        pr_before(:,(i*zs-(zs-1)):(i*zs)) = neuron.C_raw(:,z.beforepressrand{3}(i,1):z.beforepressrand{3}(i,2));
    end
    pr_before = pr_before'; % makes plotting easier

    for i = 1:size(z.afterpressrand{3},1)
            pr_after(:,(i*zs-(zs-1)):(i*zs)) = neuron.C_raw(:,z.afterpressrand{3}(i,1):z.afterpressrand{3}(i,2));
    end
    pr_after = pr_after';

    % find max, min, and increment that gives 100 thresholds between the two
    highestpointpr = max([pr_before(:);pr_after(:)]);
    lowestpointpr = min([pr_before(:);pr_after(:)]);
    incrementpr = (highestpointpr - lowestpointpr)/100; 
    
    % loop through neurons
    for j = 1:size(neuron.C_raw,1)

        combinedatapr = [pr_before(:,j);pr_after(:,j)];
        codespr = [zeros(numel(pr_before(:,j)),1);ones(numel(pr_after(:,j)),1)]; 

        %loop through thresholds
        for i = 1:100
            threshold = i*incrementpr + lowestpointpr;
            result = gt(combinedatapr,threshold);
            %calculate true positives and false positives
            tppr(i)=sum((result==1)&(codespr==1))/numel(pr_after(:,j));
            fppr(i)=sum((result==1)&(codespr==0))/numel(pr_before(:,j));
        end

        % Calculate the Riemann sum underneath the curve (AUC)
        [a,b] = sort(fppr); 
        auROCpr(f,j) = trapz(sort(fppr),tppr(b));

    end

    % sort auROC
    auROCpr_inc(f,:) = sort(auROCpr(f,:));
    

    % repeat auROC analysis on randomized head entry times
    hr_before = [];
    hr_after = [];
    for i = 1:size(z.beforeheadrand{3},1)
        hr_before(:,(i*zs-(zs-1)):(i*zs)) = neuron.C_raw(:,z.beforeheadrand{3}(i,1):z.beforeheadrand{3}(i,2));
    end
    hr_before = hr_before'; % makes plotting easier

    for i = 1:size(z.afterheadrand{3},1)
            hr_after(:,(i*zs-(zs-1)):(i*zs)) = neuron.C_raw(:,z.afterheadrand{3}(i,1):z.afterheadrand{3}(i,2));
    end
    hr_after = hr_after';

    % find max, min, and increment that gives 100 thresholds between the two
    highestpointhr = max([hr_before(:);hr_after(:)]);
    lowestpointhr = min([hr_before(:);hr_after(:)]);
    incrementhr = (highestpointhr - lowestpointhr)/100; 
    
    % loop through neurons
    for j = 1:size(neuron.C_raw,1)

        combinedatahr = [hr_before(:,j);hr_after(:,j)];
        codeshr = [zeros(numel(hr_before(:,j)),1);ones(numel(hr_after(:,j)),1)]; 

        %loop through thresholds
        for i = 1:100
            threshold = i*incrementhr + lowestpointhr;
            result = gt(combinedatahr,threshold);
            %calculate true positives and false positives
            tphr(i)=sum((result==1)&(codeshr==1))/numel(hr_after(:,j));
            fphr(i)=sum((result==1)&(codeshr==0))/numel(hr_before(:,j));
        end

        % Calculate the Riemann sum underneath the curve (AUC)
        [a,b] = sort(fphr); 
        auROChr(f,j) = trapz(sort(fphr),tphr(b));

    end

    % sort auROC
    auROChr_inc(f,:) = sort(auROChr(f,:));



% repeat auROC analysis on randomized rewarded entry times
    rr_before = [];
    rr_after = [];
    for i = 1:size(z.beforerewardrand{3},1)
        rr_before(:,(i*zs-(zs-1)):(i*zs)) = neuron.C_raw(:,z.beforerewardrand{3}(i,1):z.beforerewardrand{3}(i,2));
    end
    rr_before = rr_before'; % makes plotting easier

    for i = 1:size(z.afterrewardrand{3},1)
            rr_after(:,(i*zs-(zs-1)):(i*zs)) = neuron.C_raw(:,z.afterrewardrand{3}(i,1):z.afterrewardrand{3}(i,2));
    end
    rr_after = rr_after';

    % find max, min, and increment that gives 100 thresholds between the two
    highestpointrr = max([rr_before(:);rr_after(:)]);
    lowestpointrr = min([rr_before(:);rr_after(:)]);
    incrementrr = (highestpointrr - lowestpointrr)/100; 
    
    % loop through neurons
    for j = 1:size(neuron.C_raw,1)

        combinedatarr = [rr_before(:,j);rr_after(:,j)];
        codesrr = [zeros(numel(rr_before(:,j)),1);ones(numel(rr_after(:,j)),1)]; 

        %loop through thresholds
        for i = 1:100
            threshold = i*incrementrr + lowestpointrr;
            result = gt(combinedatarr,threshold);
            %calculate true positives and false positives
            tprr(i)=sum((result==1)&(codesrr==1))/numel(rr_after(:,j));
            fprr(i)=sum((result==1)&(codesrr==0))/numel(rr_before(:,j));
        end

        % Calculate the Riemann sum underneath the curve (AUC)
        [a,b] = sort(fprr); 
        auROCrr(f,j) = trapz(sort(fprr),tprr(b));

    end

% sort auROC
    auROCrr_inc(f,:) = sort(auROCrr(f,:));
    
% repeat auROC analysis on randomized unrewarded entry times
    ur_before = [];
    ur_after = [];
    for i = 1:size(z.beforeunrewardrand{3},1)
        ur_before(:,(i*zs-(zs-1)):(i*zs)) = neuron.C_raw(:,z.beforeunrewardrand{3}(i,1):z.beforeunrewardrand{3}(i,2));
    end
    ur_before = ur_before'; % makes plotting easier

    for i = 1:size(z.afterunrewardrand{3},1)
            ur_after(:,(i*zs-(zs-1)):(i*zs)) = neuron.C_raw(:,z.afterunrewardrand{3}(i,1):z.afterunrewardrand{3}(i,2));
    end
    ur_after = ur_after';

    % find max, min, and increment that gives 100 thresholds between the two
    highestpointur = max([ur_before(:);ur_after(:)]);
    lowestpointur = min([ur_before(:);ur_after(:)]);
    incrementur = (highestpointur - lowestpointur)/100; 
    % loop through neurons
    for j = 1:size(neuron.C_raw,1)

        combinedataur = [ur_before(:,j);ur_after(:,j)];
        codesur = [zeros(numel(ur_before(:,j)),1);ones(numel(ur_after(:,j)),1)]; 

        %loop through thresholds
        for i = 1:100
            threshold = i*incrementur + lowestpointur;
            result = gt(combinedataur,threshold);
            %calculate true positives and false positives
            tpur(i)=sum((result==1)&(codesur==1))/numel(ur_after(:,j));
            fpur(i)=sum((result==1)&(codesur==0))/numel(ur_before(:,j));
        end

        % Calculate the Riemann sum underneath the curve (AUC)
        [a,b] = sort(fpur); 
        auROCur(f,j) = trapz(sort(fpur),tpur(b));

    end

    % sort auROC
    auROCur_inc(f,:) = sort(auROCur(f,:));
end

%% find 5th and 95th percentile of randomized data
dr_5p = prctile(auROCdr_inc, 5, 1);
dr_95p = prctile(auROCdr_inc, 95, 1);

pr_5p = prctile(auROCpr_inc, 5, 1);
pr_95p = prctile(auROCpr_inc, 95, 1);

hr_5p = prctile(auROChr_inc, 5, 1);
hr_95p = prctile(auROChr_inc, 95, 1);

rr_5p = prctile(auROCrr_inc, 5, 1);
rr_95p = prctile(auROCrr_inc, 95, 1);

ur_5p = prctile(auROCur_inc, 5, 1);
ur_95p = prctile(auROCur_inc, 95, 1);


%% compare randomized to original auROC

figure; subplot(2,2,1); plot(mean(auROCpr_inc,1),'linewidth',2); 
hold on; plot(pr_5p, '--r'); hold on; plot(pr_95p,'--r');
hold on; plot(auROC_inc, '.-k');  title(sprintf('%d lever presses',size(z.rlev_press{3},1))); ylim([0,1]);

subplot(2,2,2); plot(mean(auROCdr_inc,1),'linewidth',2); 
hold on; plot(dr_5p, '--r'); hold on; plot(dr_95p,'--r');
hold on; plot(auROCd_inc, '.-k'); title(sprintf('%d dipper presentations',size(z.dip_present{3},1))); ylim([0,1]);

subplot(2,2,3); plot(mean(auROCrr_inc,1), 'linewidth',2); 
hold on; plot(rr_5p, '--r'); hold on; plot(rr_95p, '--r');
hold on; plot(auROCr_inc, '.-k'); title(sprintf('%d rewarded head entries', size(z.rewarded_entries{3},1))); ylim([0,1.0]);

subplot(2,2,4); plot(mean(auROCur_inc,1), 'linewidth',2); 
hold on; plot(ur_5p, '--r'); hold on; plot(ur_95p, '--r');
hold on; plot(auROCu_inc, '.-k'); title(sprintf('%d unrewarded head entries', size(z.unrewarded_entries{3},1))); ylim([0,1]);

set(gcf,'position',widegraph)
saveas(gcf,'randomized_v_real.png');

%% average of all zones before & after head entries

h_before1z = [];
h_before_navg = mean(h_before,2);
for i = 1:(size(h_before,1)/zs)
    h_before1z(i,:) = h_before_navg(i*zs-(zs-1):i*zs,1);
end
h_before_avg = mean(h_before1z,1);
h_after_navg = mean(h_after,2);

h_after1z = [];
for i = 1:(size(h_after,1)/zs)
    h_after1z(i,:) = h_after_navg(i*zs-(zs-1):i*zs,1);
end
h_after_avg = mean(h_after1z,1);
figure; subplot(1,3,1); 
plot([h_before_avg';h_after_avg']); title(sprintf('average zones before\n & after head entries'));
ylim([0,0.3]);
hold on; plot([zs+0.5,zs+0.5],[0,0.3], '--');

% average of all zones before & after lever presses

c_before1z = [];
c_before_navg = mean(p_before,2);
for i = 1:(size(p_before,1)/zs)
    c_before1z(i,:) = c_before_navg(i*zs-(zs-1):i*zs,1);
end
c_before_avg = mean(c_before1z,1);
c_after_navg = mean(p_after,2);

c_after1z = [];
for i = 1:(size(p_after,1)/zs)
    c_after1z(i,:) = c_after_navg(i*zs-(zs-1):i*zs,1);
end
c_after_avg = mean(c_after1z,1);
subplot(1,3,2); 
plot([c_before_avg';c_after_avg']); title(sprintf('average zones before\n & after lever presses'));
ylim([0,0.3]);
hold on; plot([zs+0.5,zs+0.5],[0,0.3], '--');
% average of all zones before & after dipper presentations

d_before1z = [];
d_before_navg = mean(d_before,2);
for i = 1:(size(d_before,1)/zs)
    d_before1z(i,:) = d_before_navg(i*zs-(zs-1):i*zs,1);
end
d_before_avg = mean(d_before1z,1);
d_after_navg = mean(d_after,2);

d_after1z = [];
for i = 1:(size(d_after,1)/zs)
    d_after1z(i,:) = d_after_navg(i*zs-(zs-1):i*zs,1);
end
d_after_avg = mean(d_after1z,1);
subplot(1,3,3); 
plot([d_before_avg';d_after_avg']); title(sprintf('average zones before\n & after dipper presentations'));
ylim([0,0.3]);
hold on; plot([zs+0.5,zs+0.5],[0,0.3], '--');
set(gcf,'position',widegraph)
saveas(gcf,'average_zones.png');

%% Activity of highest & lowest neuron in lever presses auROC

if allfigures
    figure; 
    subplot(3,1,1); plot(smoothed_rawdata(auROC==max(auROC),:)'); 
    title('highest lever press auROC neuron'); hold on; 
    for i = 1:size(z.rlev_press{3},1)
        plot([z.rlev_press{3}(i,1),z.rlev_press{3}(i,1)],[-2,7], '--k', 'linewidth', 0.8);
    end

    subplot(3,1,2);plot(smoothed_rawdata(auROC==min(auROC),:)'); 
    title('lowest lever press auROC neuron'); hold on; 
    for i = 1:size(z.rlev_press{3},1)
        plot([z.rlev_press{3}(i,1),z.rlev_press{3}(i,1)],[-2,7], '--k', 'linewidth', 0.8);
    end

    % highest neuron in lever press auROC averaged over all zones
    phi_before1z = [];
    for i = 1:(size(p_before,1)/zs)
        phi_before1z(i,:) = p_before(i*zs-(zs-1):i*zs,find(auROC==max(auROC)));
    end
    phi_before_avg = mean(phi_before1z,1);

    phi_after1z = [];
    for i = 1:(size(h_after,1)/zs)
        phi_after1z(i,:) = h_after(i*zs-(zs-1):i*zs,find(auROC==max(auROC)));
    end
    phi_after_avg = mean(phi_after1z,1);
    subplot(3,1,3); plot([phi_before_avg';phi_after_avg']); title('average activity before & after lever press of highest auROC neuron');
    hold on; plot([zs+0.5,zs+0.5],[-0.3,0.5], '--');

    plo_before1z = [];
    for i = 1:(size(p_before,1)/zs)
        plo_before1z(i,:) = p_before(i*zs-(zs-1):i*zs,find(auROC==min(auROC)));
    end
    plo_before_avg = mean(plo_before1z,1);

    plo_after1z = [];
    for i = 1:(size(h_after,1)/zs)
        plo_after1z(i,:) = h_after(i*zs-(zs-1):i*zs,find(auROC==min(auROC)));
    end
    plo_after_avg = mean(plo_after1z,1);
    subplot(3,1,3); plot([plo_before_avg';plo_after_avg']); title('average activity before & after lever press of lowest auROC neuron');
    hold on; plot([zs+0.5,zs+0.5],[-0.4,1], '--');

    %set(gcf,'position',tallgraph);
    saveas(gcf,'average_before_after_leverpress.png');
end

%% Activity of highest & lowest neuron in dipper presentation auROC
if allfigures

    figure; subplot(3,1,1); plot(smoothed_rawdata(auROCd==max(auROCd),:)'); 
    title('highest dipper presentation auROC neuron'); hold on; 
    for i = 1:size(z.dip_present{3},1)
        plot([z.dip_present{3}(i,1),z.dip_present{3}(i,1)],[-2,7], '--k', 'linewidth', 0.8);
    end

     subplot(3,1,2); plot(smoothed_rawdata(auROCd==min(auROCd),:)'); 
    title('lowest dipper presentation auROC neuron'); hold on; 
    for i = 1:size(z.dip_present{3},1)
        plot([z.dip_present{3}(i,1),z.dip_present{3}(i,1)],[-2,7], '--k', 'linewidth', 0.8);
    end

    dhi_before1z = [];
    for i = 1:(size(d_before,1)/zs)
        dhi_before1z(i,:) = d_before(i*zs-(zs-1):i*zs,find(auROCd==max(auROCd)));
    end
    dhi_before_avg = mean(dhi_before1z,1);

    dhi_after1z = [];
    for i = 1:(size(d_after,1)/zs)
        dhi_after1z(i,:) = d_after(i*zs-(zs-1):i*zs,find(auROCd==max(auROCd)));
    end
    dhi_after_avg = mean(dhi_after1z,1);
    subplot(3,1,3); plot([dhi_before_avg';dhi_after_avg']); title('average activity before & after dipper presentations of highest auROC neuron');
    hold on; plot([zs+0.5,zs+0.5],[-0.5,2.5], '--');

    dlo_before1z = [];
    for i = 1:(size(d_before,1)/zs)
        dlo_before1z(i,:) = d_before(i*zs-(zs-1):i*zs,find(auROCd==min(auROCd)));
    end
    dlo_before_avg = mean(dlo_before1z,1);

    dlo_after1z = [];
    for i = 1:(size(d_after,1)/zs)
        dlo_after1z(i,:) = d_after(i*zs-(zs-1):i*zs,find(auROCd==min(auROCd)));
    end
    dlo_after_avg = mean(dlo_after1z,1);
    subplot(3,1,3); plot([dlo_before_avg';dlo_after_avg']); title('average activity before & after dipper presentations of lowest auROC neuron');
    hold on; plot([zs+0.5,zs+0.5],[-0.6,1], '--');
    %set(gcf,'position',tallgraph);
    saveas(gcf,'average_before_after_dipper.png');

end

%% Activity of highest & lowest neuron in rewarded head entries auROC
if allfigures
    
    figure;
    subplot(3,1,1); plot(smoothed_rawdata(auROCr==max(auROCr),:)'); 
    title('highest rewarded head entry auROC neuron'); hold on; 
    for i = 1:size(z.rewarded_entries{3},1)
        plot([z.rewarded_entries{3}(i,1),z.rewarded_entries{3}(i,1)],[-2,7], '--k', 'linewidth', 0.8);
    end

    subplot(3,1,2); plot(smoothed_rawdata(auROCr==min(auROCr),:)'); 
    title('lowest rewarded head entry auROC neuron'); hold on; 
    for i = 1:size(z.rewarded_entries{3},1)
        plot([z.rewarded_entries{3}(i,1),z.rewarded_entries{3}(i,1)],[-2,7], '--k', 'linewidth', 0.8);
    end

    % highest neuron in rewarded head entries auROC averaged over all zones
    rhi_before1z = [];
    for i = 1:(size(r_before,1)/zs)
        rhi_before1z(i,:) = r_before(i*zs-(zs-1):i*zs,find(auROCr==max(auROCr)));
    end
    rhi_before_avg = mean(rhi_before1z,1);

    rhi_after1z = [];
    for i = 1:(size(r_after,1)/zs)
        rhi_after1z(i,:) = r_after(i*zs-(zs-1):i*zs,find(auROCr==max(auROCr)));
    end
    rhi_after_avg = mean(rhi_after1z,1);
    subplot(3,1,3); plot([rhi_before_avg';rhi_after_avg']); title('average activity before & after rewarded head entries of highest auROC neuron');
    hold on; plot([zs+0.5,zs+0.5],[-0.5,3.5], '--');

    % lowest neuron in rewarded head entries auROC averaged over all zones
    rlo_before1z = [];
    for i = 1:(size(r_before,1)/zs)
        rlo_before1z(i,:) = r_before(i*zs-(zs-1):i*zs,find(auROCr==min(auROCr)));
    end
    rlo_before_avg = mean(rlo_before1z,1);

    rlo_after1z = [];
    for i = 1:(size(r_after,1)/zs)
        rlo_after1z(i,:) = r_after(i*zs-(zs-1):i*zs,find(auROCr==min(auROCr)));
    end
    rlo_after_avg = mean(rlo_after1z,1);
    subplot(3,1,3); plot([rlo_before_avg';rlo_after_avg']); title('average activity before & after rewarded head entries of lowest auROC neuron');
    hold on; plot([zs+0.5,zs+0.5],[-0.6,1.4], '--');
    %set(gcf,'position',tallgraph);
    saveas(gcf,'average_before_after_rewarded_he.png');

end
%% Activity of highest & lowest neuron in unrewarded head entries auROC

if allfigures

    figure; subplot(3,1,1); plot(smoothed_rawdata(auROCu==max(auROCu),:)'); 
    title('highest unrewarded head entry auROC neuron'); hold on; 
    for i = 1:size(z.unrewarded_entries{3},1)
        plot([z.unrewarded_entries{3}(i,1),z.unrewarded_entries{3}(i,1)],[-2,7], '--k', 'linewidth', 0.8);
    end

    subplot(3,1,2); plot(smoothed_rawdata(auROCu==min(auROCu),:)'); 
    title('lowest unrewarded head entry auROC neuron'); hold on; 
    for i = 1:size(z.unrewarded_entries{3},1)
        plot([z.unrewarded_entries{3}(i,1),z.unrewarded_entries{3}(i,1)],[-2,7], '--k', 'linewidth', 0.8);
    end

    % highest neuron in unrewarded head entries auROC averaged over all zones
    uhi_before1z = [];
    for i = 1:(size(u_before,1)/zs)
        uhi_before1z(i,:) = u_before(i*zs-(zs-1):i*zs,find(auROCu==max(auROCu)));
    end
    uhi_before_avg = mean(uhi_before1z,1);

    uhi_after1z = [];
    for i = 1:(size(u_after,1)/zs)
        uhi_after1z(i,:) = u_after(i*zs-(zs-1):i*zs,find(auROCu==max(auROCu)));
    end
    uhi_after_avg = mean(uhi_after1z,1);
    subplot(3,1,3); plot([uhi_before_avg';uhi_after_avg']); title('average activity before & after unrewarded head entries of highest auROC neuron');
    hold on; plot([zs+0.5,zs+0.5],[-0.4,1], '--');

    % lowest neuron in unrewarded head entries auROC averaged over all zones
    ulo_before1z = [];
    for i = 1:(size(u_before,1)/zs)
        ulo_before1z(i,:) = u_before(i*zs-(zs-1):i*zs,find(auROCu==min(auROCu)));
    end
    ulo_before_avg = mean(ulo_before1z,1);

    ulo_after1z = [];
    for i = 1:(size(u_after,1)/zs)
        ulo_after1z(i,:) = u_after(i*zs-(zs-1):i*zs,find(auROCu==min(auROCu)));
    end
    ulo_after_avg = mean(ulo_after1z,1);
    subplot(3,1,3); plot([ulo_before_avg';ulo_after_avg']); title('average activity before & after unrewarded head entries of lowest auROC neuron');
    hold on; plot([zs+0.5,zs+0.5],[-0.2,1.4], '--');
    %set(gcf,'position',tallgraph);
    saveas(gcf,'average_before_after_unrewarded_he.png');

end

%% find location of each neuron
neuron_image=[];

for i = 1:size(neuron.A,2)
    neuron_image(:,:,i) = reshape(neuron.A(:,i),[size(neuron.A,1)^(1/2),size(neuron.A,1)^(1/2)]);
end

n_image = double(sum(neuron_image,3)>0);

%% locate top 10 and bottom 10 neurons of each auROC

% top10 anonymous
top10 = @(index) index(size(index,2)-10:size(index,2));
% onesigmaabove anonymous (returns indices of aurocs one sigma above their
% distribution)
onesigmaabove = @(aurocs) find(aurocs>(mean(aurocs)+std(aurocs)))

% bottom10 anonymous
bottom10 = @(index) index(1:10);
% onesigmabelow anonymous (returns indices of aurocs one sigma below their
% distribution)
onesigmabelow = @(aurocs) find(aurocs<(mean(aurocs)-std(aurocs)))

% showclusters anonymous
showclusters = @(neuron_image,indexhigh,indexlow) max(neuron_image(:,:,indexhigh)>0,[],3)-max(neuron_image(:,:,indexlow)>0,[],3);

% Calculate top and bottom groups for all behaviors
% lever presses
top10_lp = top10(index); bottom10_lp = bottom10(index);
top10_r = top10(indexr); bottom10_r = bottom10(indexr);
top10_d = top10(indexd); bottom10_d = bottom10(indexd);
top10_u = top10(indexu); bottom10_u = bottom10(indexu);

% lever presses
figure; 
subplot(2,2,1); imagesc(showclusters(neuron_image,top10_lp,bottom10_lp)); title(sprintf('Top 10 for day %1.0f and time = %1.0f s\n \n Lever press neurons',day,zs/5));
% dipper presentations
subplot(2,2,2); imagesc(showclusters(neuron_image,top10_d,bottom10_d)); title(sprintf('[BLUE=inhibited; ORANGE=excited]\n \n Dipper presentation neurons'));
% rewarded entries
subplot(2,2,3); imagesc(showclusters(neuron_image,top10_r,bottom10_r)); title('Rewarded entry neurons');
% unrewarded entries
subplot(2,2,4); imagesc(showclusters(neuron_image,top10_u,bottom10_u)); title('Unrewarded entry neurons');

for i = 1:4; subplot(2,2,i); colormap(jet(5)); end

%%
% Outputs to SAVE to a master structure (called data_by_day)

data_by_day(day,zs).auROC = auROC;
data_by_day(day,zs).auROCd = auROCd;
data_by_day(day,zs).auROCu = auROCu;
data_by_day(day,zs).auROCr = auROCr;

data_by_day(day,zs).dr_5p=dr_5p;
data_by_day(day,zs).dr_95p=dr_95p;
data_by_day(day,zs).pr_5p=pr_5p;
data_by_day(day,zs).pr_95p=pr_95p;
data_by_day(day,zs).hr_5p=hr_5p;
data_by_day(day,zs).hr_95p=hr_95p;
data_by_day(day,zs).rr_5p=rr_5p;
data_by_day(day,zs).rr_95p=rr_95p;
data_by_day(day,zs).ur_5p=ur_5p;
data_by_day(day,zs).ur_95p=ur_95p;

data_by_day(day,zs).auROCpr_inc = auROCpr_inc;
data_by_day(day,zs).auROCdr_inc = auROCdr_inc;
data_by_day(day,zs).auROChr_inc = auROChr_inc;
data_by_day(day,zs).auROCrr_inc = auROCrr_inc;
data_by_day(day,zs).auROCur_inc = auROCur_inc;

data_by_day(day,zs).cluster_lp = showclusters(neuron_image,top10_lp,bottom10_lp);
data_by_day(day,zs).cluster_d = showclusters(neuron_image,top10_lp,bottom10_d);
data_by_day(day,zs).cluster_r = showclusters(neuron_image,top10_lp,bottom10_r);
data_by_day(day,zs).cluster_u = showclusters(neuron_image,top10_lp,bottom10_u);

data_by_day(day,zs).corrcoef = corrcoef([auROC(indexd);auROCd(indexd);auROCr(indexd);auROCu(indexd)]');

