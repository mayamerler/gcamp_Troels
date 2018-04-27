%% Loads data for mouse 7677 over 5 days

files{1} = 'D:\troels gcamp\Workspaces\Workspace 7677 RR5 Day1 calculated.mat';
files{2} = 'D:\troels gcamp\Workspaces\Workspace 7677 RR5 Day2 calculated.mat';
files{3} = 'D:\troels gcamp\Workspaces\Workspace 7677 RR5 Day3 calculated.mat';
files{4} = 'D:\troels gcamp\Workspaces\Workspace 7677 RR5 Day4 calculated.mat';
files{5} = 'D:\troels gcamp\Workspaces\Workspace 7677 RR5 Day5 original.mat';

% List of figures generated:

% Let's fill this in...
%
%

% Other needed edits:
%
% (3)  
%
% Use allfigures = 0 to suppress figures (only analyze data)

% List of dependencies:
%
% Other functions needed for this to run, e.g., mayabeforezones (which
% themselves need to be a bit more commented)


% Parameters for analysis
day=5; % Set day to a value between 1 and 5, run the whole script to generate figures

load(files{day});
load('D:\troels gcamp\Workspaces\rr5data.mat'); % This is to get the RR5 data without Excel 

zs=15; % sets size of time zone in samples (sampling rate for this data is 5 Hz

widegraph = [45 70 1000 550];

%% Initializing 

set(0,'defaultlinelinewidth',1.5);
set(0,'defaultfigurecolor','w');

allfigures=0;
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
rr5_7677 = rr5data(day).m7677;
rr5_7677(find(rr5_7677(:,1)==0),:)=[]; % Eliminate random zeroes

rr5head = rr5_7677(:,1)*5; %converts seconds to frames and caps at 1500 frames
rr5dipper = rr5_7677(:,2)*5;
rr5press = rr5_7677(:,3)*5;
z.head_entries = rr5head(find(rr5head<1500));
z.dip_present = rr5dipper(find(rr5dipper<1500));
z.rlev_press = rr5press(find(rr5press<1500));

%% separate rewarded from non-rewarded head entries

rcounter = 1;
ucounter = 1;
for i = 1:size(z.head_entries,1)
    if sum((z.head_entries(i) > z.dip_present) & (z.head_entries(i) < (z.dip_present+25)))>0
        z.rewarded_entries(rcounter,1) = z.head_entries(i);
        rcounter = rcounter + 1;
    else
        z.unrewarded_entries(ucounter,1) = z.head_entries(i);
        ucounter = ucounter + 1;
    end
end

%% finding time zones before and after each event

z.beforepress = mayabeforezones(z.rlev_press,zs);
z.afterpress = mayaafterzones(z.rlev_press,zs);
z.beforedip = mayabeforezones(z.dip_present,zs);
z.afterdip = mayaafterzones(z.dip_present,zs);
z.beforehead = mayabeforezones(z.head_entries,zs);
z.afterhead = mayaafterzones(z.head_entries,zs);
z.beforereward = mayabeforezones(z.rewarded_entries,zs);
z.afterreward = mayaafterzones(z.rewarded_entries,zs);
z.beforeunreward = mayabeforezones(z.unrewarded_entries,zs);
z.afterunreward = mayaafterzones(z.unrewarded_entries,zs);

%% auROC per event (rewarded or unrewarded entries)
for i = 1:size(z.afterreward,1)
    auROCr_event(i,:) = mayaauroc(neuron.C_raw,z.beforereward(i,1):z.beforereward(i,2),z.afterreward(i,1):z.afterreward(i,2));
end
[auROCr_event_inc, indexr_event] = sort(auROCr_event'); 
for i = 1:size(z.afterunreward,1)
    auROCu_event(i,:) = mayaauroc(neuron.C_raw,z.beforeunreward(i,1):z.beforeunreward(i,2),z.afterunreward(i,1):z.afterunreward(i,2));
end
[auROCu_event_inc, indexu_event] = sort(auROCu_event');

%% rank each cell over days
%rewarded
clearvars cell_rank* rank_mean* rank_std*
for i = 1:size(auROCr_event,2) 
    [cell_rankr(:,i),b] = find(indexr_event==i);
end

for i = 1:size(cell_rankr,2)
   rank_meanr(i) = mean(cell_rankr(:,i));
   rank_stdr(i) = std(cell_rankr(:,i));
end
[rank_meanr_inc,mean_indexr] = sort(rank_meanr);
[fr,xr] = ecdf(rank_stdr);
centerr = mean(rank_meanr);
dist_from_centerr = abs(rank_meanr - centerr);
correlationr = corrcoef(dist_from_centerr,rank_stdr);
fitr = polyfit(dist_from_centerr,rank_stdr,1);

%unrewarded
for i = 1:size(auROCu_event,2) 
    [cell_ranku(:,i),b] = find(indexu_event==i);
end

for i = 1:size(cell_ranku,2)
   rank_meanu(i) = mean(cell_ranku(:,i));
   rank_stdu(i) = std(cell_ranku(:,i));
end
[rank_meanu_inc,mean_indexu] = sort(rank_meanu);
[fu,xu] = ecdf(rank_stdu);
centeru = mean(rank_meanu);
dist_from_centeru = abs(rank_meanu - centeru);
correlationu = corrcoef(dist_from_centeru,rank_stdu);
fitu = polyfit(dist_from_centeru,rank_stdu,1);

if allfigures; figure; 
subplot(2,2,1); imagesc(cell_rankr(:,mean_indexr)'); title('rewarded'); xlabel('event'); ylabel('neuron');
subplot(2,2,2); imagesc(cell_ranku(:,mean_indexu)'); title('unrewarded'); xlabel('event'); ylabel('neuron');
subplot(2,2,3); plot(dist_from_centerr,rank_stdr,'o'); hold on; 
plot([0:max(dist_from_centerr+1)],([0:max(dist_from_centerr+1)]*fitr(1)+fitr(2))); title(sprintf('correlation coefficient = %1.3f',correlationr(2)));
xlabel('distance of cell mean from overall mean'); ylabel('standard deviation'); 
xlim([0,centerr]); ylim([0,max([rank_stdr,rank_stdu])+1]);
subplot(2,2,4); plot(dist_from_centeru,rank_stdu,'o'); hold on; 
plot([0:max(dist_from_centeru+1)],([0:max(dist_from_centeru+1)]*fitu(1)+fitu(2))); title(sprintf('correlation coefficient = %1.3f',correlationu(2)));
xlabel('distance of cell mean from overall mean'); ylabel('standard deviation'); 
xlim([0,centerr]); ylim([0,max([rank_stdr,rank_stdu])+1]);
set(gcf,'position',widegraph);
name = sprintf('D:/troels gcamp/cell rank figs/7675day%1.0f.svg',day);
saveas(gcf,name);
end

%% find location of each neuron
neuron_image=[];

for i = 1:size(neuron.A,2)
    neuron_image(:,:,i) = reshape(neuron.A(:,i),[size(neuron.A,1)^(1/2),size(neuron.A,1)^(1/2)]);
end

n_image = double(sum(neuron_image,3)>0);

%% location of neurons ranked by mean rewarded entry
rank_imager = zeros(315,315,size(neuron_image,3));
for i = 1:size(rank_meanr,2)
    rank_imager(:,:,i) = ((neuron_image(:,:,i)>0).*(max(rank_meanr)+1-rank_meanr(i)));
end
if allfigures
    figure; imagesc(max(rank_imager,[],3),[0 max((max(rank_meanr)+1-rank_meanr))]); 
end
%% by unrewarded entry for comparison
rank_imageu = zeros(315,315,size(neuron_image,3));
for i = 1:size(rank_meanu,2)
    rank_imageu(:,:,i) = ((neuron_image(:,:,i)>0).*(max(rank_meanu)+1-rank_meanu(i)));
end
if allfigures
figure; imagesc(max(rank_imageu,[],3),[0 max((max(rank_meanr)+1-rank_meanr))]); 
end
%[0 max((max(rank_meanr)+1-rank_meanr))]
%% auROC for lever presses
% grouping all data from 'before' zones and all data from 'after' zones
p_before = [];
for i = 1:size(z.beforepress,1)
    p_before(:,(i*zs-(zs-1)):(i*zs)) = neuron.C_raw(:,z.beforepress(i,1):z.beforepress(i,2));
    p_index_before((i*zs-(zs-1)):(i*zs)) = z.beforepress(i,1):z.beforepress(i,2);
end
p_before = p_before'; % makes plotting easier

p_after = [];
for i = 1:size(z.afterpress,1)
    p_after(:,(i*zs-(zs-1)):(i*zs)) = neuron.C_raw(:,z.afterpress(i,1):z.afterpress(i,2));
    p_index_after((i*zs-(zs-1)):(i*zs)) = z.afterpress(i,1):z.afterpress(i,2);
end
p_after = p_after';

auROCp = [];
auROCp = mayaauroc(neuron.C_raw,p_index_before,p_index_after);

[auROCp_inc, indexp] = sort(auROCp);

%% repeat for dipper presentations
d_before = [];
for i = 1:size(z.beforedip,1)
    d_before(:,(i*zs-(zs-1)):(i*zs)) = neuron.C_raw(:,z.beforedip(i,1):z.beforedip(i,2));
    d_index_before((i*zs-(zs-1)):(i*zs)) = z.beforedip(i,1):z.beforedip(i,2);
end
d_before = d_before'; % makes plotting easier
d_after = [];
for i = 1:size(z.afterdip,1)
    d_after(:,(i*zs-(zs-1)):(i*zs)) = neuron.C_raw(:,z.afterdip(i,1):z.afterdip(i,2));
    d_index_after((i*zs-(zs-1)):(i*zs)) = z.afterdip(i,1):z.afterdip(i,2);
end
d_after = d_after';

auROCd = [];
auROCd = mayaauroc(neuron.C_raw,d_index_before,d_index_after);

[auROCd_inc, indexd] = sort(auROCd);

%% repeat for head entries

h_before = [];
h_after = [];
for i = 1:size(z.beforehead,1)
    h_before(:,(i*zs-(zs-1)):(i*zs)) = neuron.C_raw(:,z.beforehead(i,1):z.beforehead(i,2));
    h_index_before((i*zs-(zs-1)):(i*zs)) = z.beforehead(i,1):z.beforehead(i,2);
end
h_before = h_before'; % makes plotting easier

for i = 1:size(z.afterhead,1)
    h_after(:,(i*zs-(zs-1)):(i*zs)) = neuron.C_raw(:,z.afterhead(i,1):z.afterhead(i,2));
    h_index_after((i*zs-(zs-1)):(i*zs)) = z.afterhead(i,1):z.afterhead(i,2);
end
h_after = h_after';

auROCh = [];
auROCh = mayaauroc(neuron.C_raw,h_index_before,h_index_after);
[auROCh_inc,indexh] = sort(auROCh);

%% auROC for just rewarded head entries

r_before=[];
for i = 1:size(z.beforereward,1)
    r_before(:,(i*zs-(zs-1)):(i*zs)) = neuron.C_raw(:,z.beforereward(i,1):z.beforereward(i,2));
    r_index_before((i*zs-(zs-1)):(i*zs)) = z.beforereward(i,1):z.beforereward(i,2);
end
r_before = r_before'; % makes plotting easier

r_after=[];
for i = 1:size(z.afterreward,1)
    r_after(:,(i*zs-(zs-1)):(i*zs)) = neuron.C_raw(:,z.afterreward(i,1):z.afterreward(i,2));
    r_index_after((i*zs-(zs-1)):(i*zs)) = z.afterreward(i,1):z.afterreward(i,2);
end
r_after = r_after';

auROCr = [];
auROCr = mayaauroc(neuron.C_raw,r_index_before,r_index_after);
[auROCr_inc, indexr] = sort(auROCr);

%% auROC for just unrewarded head entries

u_before = [];
for i = 1:size(z.beforeunreward,1)
    u_before(:,(i*zs-(zs-1)):(i*zs)) = neuron.C_raw(:,z.beforeunreward(i,1):z.beforeunreward(i,2));
    u_index_before((i*zs-(zs-1)):(i*zs)) = z.beforeunreward(i,1):z.beforeunreward(i,2);
end
u_before = u_before'; % makes plotting easier

u_after = [];
for i = 1:size(z.afterunreward,1)
    u_after(:,(i*zs-(zs-1)):(i*zs)) = neuron.C_raw(:,z.afterunreward(i,1):z.afterunreward(i,2));
    u_index_after((i*zs-(zs-1)):(i*zs)) = z.afterunreward(i,1):z.afterunreward(i,2);
end
u_after = u_after';

auROCu = [];
auROCu = mayaauroc(neuron.C_raw,u_index_before,u_index_after);
[auROCu_inc, indexu] = sort(auROCu);
auROCu_smooth = smooth(auROCu_inc,5);

if allfigures; 
figure; plot(auROCu_inc, '.-','markersize',12);  title('auROC unrewarded/rewarded head entries');
hold on; plot(auROCr_inc, 'r.-','markersize',12);
end

%% plot auROC for all events

ylims = [0,1];

if allfigures;
figure; subplot(2,2,1); plot(auROCp_inc, 'o-'); title('lever presses'); ylim(ylims);
subplot(2,2,2); plot(auROCd_inc, 'o-');  title('dipper presentations'); ylim(ylims);
subplot(2,2,3); plot(auROCh_inc, 'o-'); title('all head entries'); ylim(ylims); 
subplot(2,2,4); plot(auROCu_inc, 'ro-'); title('head entries by reward'); ylim(ylims);hold on;
subplot(2,2,4); plot(auROCr_inc, 'go-'); ylim(ylims); 
%set(gcf,'position',boxgraph)
saveas(gcf,'auROC_by_event.png'); end

%% averaged data with stimulus info

if allfigures;
figure; subplot(3,1,1); plot(smooth_avgC, 'k'); title('averaged data with head entries'); hold on; 

for h = 1:size(z.head_entries,1)
    plot([z.head_entries(h,1),z.head_entries(h,1)],[-0.2,0.7], '--r','linewidth',0.8);
end

subplot(3,1,2); plot(smooth_avgC, 'k'); title('averaged data with lever presses'); hold on;
for p = 1:size(z.rlev_press,1)
    plot([z.rlev_press(p,1),z.rlev_press(p,1)],[-0.2,0.7],'--g','linewidth',0.8);
end

subplot(3,1,3); plot(smooth_avgC, 'k'); title('averaged data with dipper presentations'); hold on;
for i = 1:size(z.dip_present,1)
    plot([z.dip_present(i,1),z.dip_present(i,1)],[-0.2,0.7], '--b','linewidth',0.8);
end 
set(gcf,'position',widegraph)
saveas(gcf,'event_trigged_average.png');
end
%% averaged data with head entries & dipper presentations

if allfigures;
figure; plot(smooth_avgC, 'k'); title('averaged data with head entries (r) & dipper presentations (b)'); hold on;
for h = 1:size(z.head_entries,1)
    plot([z.head_entries(h,1),z.head_entries(h,1)],[-0.2,0.7], '--r','linewidth',0.8);
end
hold on;
for i = 1:size(z.dip_present,1)
    plot([z.dip_present(i,1),z.dip_present(i,1)],[-0.2,0.7], '--b','linewidth',0.8);
end 
set(gcf,'position',widegraph)
saveas(gcf,'dipper_trigged_average.png');
end

%% maintaining same order of neurons over all plots (indexed by dipper presentation order)
if allfigures;
figure; subplot(1,4,2); plot(auROCp(indexd), 'o'); title('lever presses'); ylim([0.1,0.9]);
subplot(1,4,1); plot(auROCd(indexd), 'ro');  title('dipper presentations'); ylim([0.1,0.9]);
ylabel(sprintf('Day %1.0f auROC per neuron\n (in order of increasing value for reward)\n with before and after comparison window of %1.0f s',day,zs/5))
subplot(1,4,3); plot(auROCr(indexd), 'o'); title('rewarded head entries'); ylim([0.1,0.9]);
subplot(1,4,4); plot(auROCu(indexd), 'o'); title('unrewarded head entries'); ylim([0.1,0.9]);

set(gcf,'position',widegraph); saveas(gcf,'auROC_indexed_to_dippervalues.png');
end

%% maintaining same order of neurons over all plots (indexed by rewarded entry order)
if allfigures;
figure; subplot(1,4,2); plot(auROCp(indexr), 'o'); title('lever presses'); ylim([0.1,0.9]);
subplot(1,4,3); plot(auROCd(indexr), 'o');  title('dipper presentations'); ylim([0.1,0.9]);
subplot(1,4,1); plot(auROCr(indexr), 'ro'); title('rewarded head entries'); ylim([0.1,0.9]);
ylabel(sprintf('Day %1.0f auROC per neuron\n (in order of increasing value for reward)\n with before and after comparison window of %1.0f s',day,zs/5))
subplot(1,4,4); plot(auROCu(indexr), 'o'); title('unrewarded head entries'); ylim([0.1,0.9]);

set(gcf,'position',widegraph); saveas(gcf,'auROC_indexed_to_reward.png');
end

%% maintaining same order of neurons over all plots (indexed by lever press order)
if allfigures;
figure; subplot(1,4,1); plot(auROCp(indexp), 'ro'); title('lever presses'); ylim([0.1,0.9]);
ylabel(sprintf('Day %1.0f auROC per neuron\n (in order of increasing value for reward)\n with before and after comparison window of %1.0f s',day,zs/5))
subplot(1,4,2); plot(auROCd(indexp), 'o');  title('dipper presentations'); ylim([0.1,0.9]);
subplot(1,4,3); plot(auROCr(indexp), 'o'); title('rewarded head entries'); ylim([0.1,0.9]);
subplot(1,4,4); plot(auROCu(indexp), 'o'); title('unrewarded head entries'); ylim([0.1,0.9]);
set(gcf,'position',widegraph); saveas(gcf,'auROC_indexed_to_leverpress.png');
end

%% maintaining same order of neurons over all plots (indexed by unrewarded entry order)
if allfigures;
figure; subplot(1,4,2); plot(auROCp(indexu), 'o'); title('lever presses'); ylim([0.1,0.9]);
subplot(1,4,4); plot(auROCd(indexu), 'o');  title('dipper presentations'); ylim([0.1,0.9]);
subplot(1,4,3); plot(auROCr(indexu), 'o'); title('rewarded head entries'); ylim([0.1,0.9]);
subplot(1,4,1); plot(auROCu(indexu), 'ro'); title('unrewarded head entries'); ylim([0.1,0.9]);
ylabel(sprintf('Day %1.0f auROC per neuron\n (in order of increasing value for reward)\n with before and after comparison window of %1.0f s',day,zs/5))
set(gcf,'position',widegraph); saveas(gcf,'auROC_indexed_to_reward.png');
end

%% correlation coefficients between groups
if allfigures;
figure; imagesc(corrcoef([auROCp(indexd);auROCd(indexd);auROCr(indexd);auROCu(indexd)]'));
xticks([1:4]); yticks([1:4]); xticklabels({'lever','dipper','rewarded','unrewarded'});
yticklabels({'lever','dipper','rewarded','unrewarded'}); 
title(sprintf('Correlation coefficients (neuron auROC x neuron auROC)\nfor different periods of day %1.0f\n with before/after comparison window of %1.0f s',day,zs/5));
colorbar(); end

%% Activity of two neurons with auROC over 0.7 in dipper presentation data
for j = 1:size(neuron.C_raw,1)
    smoothed_rawdata(j,:) = smooth(neuron.C_raw(j,:),5);
end

cutoff = mean(auROCd)+std(auROCd);
if allfigures; figure; subplot(2,2,1); plot(smoothed_rawdata(auROCd>cutoff,:)'); hold on; 
for i = 1:size(z.dip_present,1)
    plot([z.dip_present(i,1),z.dip_present(i,1)],[-2,7], 'k');
end
title(sprintf('%1.0f neurons with auROC>%1.2f (locked to dipper)',sum(auROCd>cutoff),cutoff));

cutoff = mean(auROCd)-std(auROCd);
subplot(2,2,2); plot(smoothed_rawdata(auROCd<cutoff,:)'); hold on; 
for i = 1:size(z.dip_present,1)
    plot([z.dip_present(i,1),z.dip_present(i,1)],[-2,7], 'k');
end
title(sprintf('%1.0f neurons with auROC<%1.2f (locked to dipper)',sum(auROCd<cutoff),cutoff));

cutoff = mean(auROCd)+std(auROCd);
subplot(2,2,3); plot(mean(smoothed_rawdata(auROCd>cutoff,:)',2)); hold on; 
for i = 1:size(z.dip_present,1)
    plot([z.dip_present(i,1),z.dip_present(i,1)],[-2,7], 'k');
end
ylim([-1,3]);
title(sprintf('%1.0f neurons with auROC>%1.2f (locked to dipper)',sum(auROCd>cutoff),cutoff));

cutoff = mean(auROCd)-std(auROCd);
subplot(2,2,4); plot(mean(smoothed_rawdata(auROCd<cutoff,:)',2)); hold on; 
for i = 1:size(z.dip_present,1)
    plot([z.dip_present(i,1),z.dip_present(i,1)],[-2,7], 'k');
end
ylim([-1,3]);
title(sprintf('%1.0f neurons with auROC<%1.2f (locked to dipper)',sum(auROCd<cutoff),cutoff));
set(gcf,'position',widegraph); saveas(gcf,'averages_by_auROC.png');
end

%% randomizing times of events 
% clearvars auROCdr auROChr auROCpr auROCrr auROCur;

for f = 1:200 % Iterations for the randomization
    
    rand_dip = sort(randi(1500, size(z.dip_present,1),1));
    rand_press = sort(randi(1500, size(z.rlev_press,1),1));
    rand_head = sort(randi(1500, size(z.head_entries,1),1));
    rand_reward = sort(randi(1500, size(z.rewarded_entries,1),1));
    rand_unreward = sort(randi(1500, size(z.unrewarded_entries,1),1));
    
% finding time zones before and after each randomized event
    
    z.beforepressrand = mayabeforezones(rand_press,zs);
    z.afterpressrand = mayaafterzones(rand_press,zs);
    
    z.beforediprand = mayabeforezones(rand_dip,zs);
    z.afterdiprand = mayaafterzones(rand_dip,zs);
    
    z.beforeheadrand = mayabeforezones(rand_head,zs);
    z.afterheadrand = mayaafterzones(rand_head,zs);
    
    z.beforerewardrand = mayabeforezones(rand_reward,zs);
    z.afterrewardrand = mayaafterzones(rand_reward,zs);
    
    z.beforeunrewardrand = mayabeforezones(rand_unreward,zs);
    z.afterunrewardrand = mayaafterzones(rand_unreward,zs);    

% repeat auROC analysis on randomized dipper presentation times
    dr_before = [];
    dr_after = [];
    for i = 1:size(z.beforediprand,1)
        dr_before(:,(i*zs-(zs-1)):(i*zs)) = neuron.C_raw(:,z.beforediprand(i,1):z.beforediprand(i,2));
        dr_index_before((i*zs-(zs-1)):(i*zs)) = z.beforediprand(i,1):z.beforediprand(i,2);
    end
    dr_before = dr_before'; % makes plotting easier

    for i = 1:size(z.afterdiprand,1)
        dr_after(:,(i*zs-(zs-1)):(i*zs)) = neuron.C_raw(:,z.afterdiprand(i,1):z.afterdiprand(i,2));
        dr_index_after((i*zs-(zs-1)):(i*zs)) = z.afterdiprand(i,1):z.afterdiprand(i,2);
    end
    dr_after = dr_after';
    
    auROCdr(f,:) = mayaauroc(neuron.C_raw,dr_index_before,dr_index_after);
    auROCdr_inc(f,:) = sort(auROCdr(f,:));
   
    
% repeat auROC analysis on randomized lever press times
    pr_before = [];
    pr_after = [];
    for i = 1:size(z.beforepressrand,1)
        pr_before(:,(i*zs-(zs-1)):(i*zs)) = neuron.C_raw(:,z.beforepressrand(i,1):z.beforepressrand(i,2));
        pr_index_before((i*zs-(zs-1)):(i*zs)) = z.beforepressrand(i,1):z.beforepressrand(i,2);
    end
    pr_before = pr_before'; % makes plotting easier

    for i = 1:size(z.afterpressrand,1)
        pr_after(:,(i*zs-(zs-1)):(i*zs)) = neuron.C_raw(:,z.afterpressrand(i,1):z.afterpressrand(i,2));
        pr_index_after((i*zs-(zs-1)):(i*zs)) = z.afterpressrand(i,1):z.afterpressrand(i,2);
    end
    pr_after = pr_after';

    auROCpr(f,:) = mayaauroc(neuron.C_raw,pr_index_before,pr_index_after);
    auROCpr_inc(f,:) = sort(auROCpr(f,:));
    

    % repeat auROC analysis on randomized head entry times
    hr_before = [];
    hr_after = [];
    for i = 1:size(z.beforeheadrand,1)
        hr_before(:,(i*zs-(zs-1)):(i*zs)) = neuron.C_raw(:,z.beforeheadrand(i,1):z.beforeheadrand(i,2));
        hr_index_before((i*zs-(zs-1)):(i*zs)) = z.beforeheadrand(i,1):z.beforeheadrand(i,2);
    end
    hr_before = hr_before'; % makes plotting easier

    for i = 1:size(z.afterheadrand,1)
        hr_after(:,(i*zs-(zs-1)):(i*zs)) = neuron.C_raw(:,z.afterheadrand(i,1):z.afterheadrand(i,2));
        hr_index_after((i*zs-(zs-1)):(i*zs)) = z.afterheadrand(i,1):z.afterheadrand(i,2);
    end
    hr_after = hr_after';

    auROChr(f,:) = mayaauroc(neuron.C_raw,hr_index_before,hr_index_after);
    auROChr_inc(f,:) = sort(auROChr(f,:));



% repeat auROC analysis on randomized rewarded entry times
    rr_before = [];
    rr_after = [];
    for i = 1:size(z.beforerewardrand,1)
        rr_before(:,(i*zs-(zs-1)):(i*zs)) = neuron.C_raw(:,z.beforerewardrand(i,1):z.beforerewardrand(i,2));
        rr_index_before((i*zs-(zs-1)):(i*zs)) = z.beforerewardrand(i,1):z.beforerewardrand(i,2);
    end
    rr_before = rr_before'; % makes plotting easier

    for i = 1:size(z.afterrewardrand,1)
        rr_after(:,(i*zs-(zs-1)):(i*zs)) = neuron.C_raw(:,z.afterrewardrand(i,1):z.afterrewardrand(i,2));
        rr_index_after((i*zs-(zs-1)):(i*zs)) = z.afterrewardrand(i,1):z.afterrewardrand(i,2);
    end
    rr_after = rr_after';
    
    auROCrr(f,:) = mayaauroc(neuron.C_raw,rr_index_before,rr_index_after);
    auROCrr_inc(f,:) = sort(auROCrr(f,:));
    
% repeat auROC analysis on randomized unrewarded entry times
    ur_before = [];
    ur_after = [];
    for i = 1:size(z.beforeunrewardrand,1)
        ur_before(:,(i*zs-(zs-1)):(i*zs)) = neuron.C_raw(:,z.beforeunrewardrand(i,1):z.beforeunrewardrand(i,2));
        ur_index_before((i*zs-(zs-1)):(i*zs)) = z.beforeunrewardrand(i,1):z.beforeunrewardrand(i,2);
    end
    ur_before = ur_before'; % makes plotting easier
    
    for i = 1:size(z.afterrewardrand,1)
        ur_after(:,(i*zs-(zs-1)):(i*zs)) = neuron.C_raw(:,z.afterunrewardrand(i,1):z.afterunrewardrand(i,2));
        ur_index_after((i*zs-(zs-1)):(i*zs)) = z.afterunrewardrand(i,1):z.afterunrewardrand(i,2);
    end
    ur_after = ur_after';
    
    auROCur(f,:) = mayaauroc(neuron.C_raw,ur_index_before,ur_index_after);
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
if allfigures;
figure; subplot(2,2,1); plot(mean(auROCpr_inc,1),'linewidth',2); 
hold on; plot(pr_5p, '--r'); hold on; plot(pr_95p,'--r');
hold on; plot(auROCp_inc, '.-k');  title(sprintf('%d lever presses',size(z.rlev_press,1))); ylim([0,1]);

subplot(2,2,2); plot(mean(auROCdr_inc,1),'linewidth',2); 
hold on; plot(dr_5p, '--r'); hold on; plot(dr_95p,'--r');
hold on; plot(auROCd_inc, '.-k'); title(sprintf('%d dipper presentations',size(z.dip_present,1))); ylim([0,1]);

subplot(2,2,3); plot(mean(auROCrr_inc,1), 'linewidth',2); 
hold on; plot(rr_5p, '--r'); hold on; plot(rr_95p, '--r');
hold on; plot(auROCr_inc, '.-k'); title(sprintf('%d rewarded head entries', size(z.rewarded_entries,1))); ylim([0,1.0]);

subplot(2,2,4); plot(mean(auROCur_inc,1), 'linewidth',2); 
hold on; plot(ur_5p, '--r'); hold on; plot(ur_95p, '--r');
hold on; plot(auROCu_inc, '.-k'); title(sprintf('%d unrewarded head entries', size(z.unrewarded_entries,1))); ylim([0,1]);

set(gcf,'position',widegraph); saveas(gcf,'randomized_v_real.png');
end

%% Percent of points outside confidence interval
count = 0;
for i = 1:size(auROCr_inc,2)
    if auROCr_inc(i) < rr_5p(i)
        count = count+1;
        ci_distr(count) = rr_5p(i)-auROCr_inc(i);
    end
    if auROCr_inc(i) > rr_95p(i)
        count = count + 1;
        ci_distr(count) = auROCr_inc(i) - rr_5p(i);
    end
end
ci_outsider = count/i;
net_distr = sum(ci_distr);

count = 0;
for i = 1:size(auROCu_inc,2)
    if auROCu_inc(i) < ur_5p(i)
        count = count + 1;        
        ci_distu(count) = ur_5p(i)-auROCu_inc(i);
    end
    if auROCu_inc(i) > ur_95p(i)
        count = count + 1;        
        ci_distu(count) = auROCu_inc(i) - ur_5p(i);
    end
end
ci_outsideu = count/i;
net_distu = sum(ci_distu);

count = 1;
for i = 1:size(auROCu_inc,2)
    if auROCu_inc(i) < ur_5p(i)
        ci_distu(count) = ur_5p(i)-auROCu_inc(i);
        count = count+1;
    end
    if auROCu_inc(i) > ur_95p(i)
        ci_distu(count) = auROCu_inc(i) - ur_5p(i);
        count = count + 1;
    end
end
ci_outsideu = count/i;
net_distu = sum(ci_distu);

count = 1;
for i = 1:size(auROCu_inc,2)
    if auROCu_inc(i) < ur_5p(i)
        ci_distu(count) = ur_5p(i)-auROCu_inc(i);
        count = count+1;
    end
    if auROCu_inc(i) > ur_95p(i)
        ci_distu(count) = auROCu_inc(i) - ur_5p(i);
        count = count + 1;
    end
end
ci_outsideu = count/i;
net_distu = sum(ci_distu);

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

if allfigures; figure; subplot(1,3,1); plot([h_before_avg';h_after_avg']); 
title(sprintf('average zones before\n & after head entries')); ylim([0,0.3]);
hold on; plot([zs+0.5,zs+0.5],[0,0.3], '--');
subplot(1,3,2); plot([c_before_avg';c_after_avg']); 
title(sprintf('average zones before\n & after lever presses')); ylim([0,0.3]);
hold on; plot([zs+0.5,zs+0.5],[0,0.3], '--');
subplot(1,3,3); plot([d_before_avg';d_after_avg']); 
title(sprintf('average zones before\n & after dipper presentations')); ylim([0,0.3]);
hold on; plot([zs+0.5,zs+0.5],[0,0.3], '--');
set(gcf,'position',widegraph); saveas(gcf,'average_zones.png');
end

%% Activity of highest & lowest neuron in lever presses auROC

if allfigures
    figure; 
    subplot(3,1,1); plot(smoothed_rawdata(auROCp==max(auROCp),:)'); 
    title('highest lever press auROC neuron'); hold on; 
    for i = 1:size(z.rlev_press,1)
        plot([z.rlev_press(i,1),z.rlev_press(i,1)],[-2,7], '--k', 'linewidth', 0.8);
    end

    subplot(3,1,2);plot(smoothed_rawdata(auROCp==min(auROCp),:)'); 
    title('lowest lever press auROC neuron'); hold on; 
    for i = 1:size(z.rlev_press,1)
        plot([z.rlev_press(i,1),z.rlev_press(i,1)],[-2,7], '--k', 'linewidth', 0.8);
    end

    % highest neuron in lever press auROC averaged over all zones
    phi_before1z = [];
    for i = 1:(size(p_before,1)/zs)
        phi_before1z(i,:) = p_before(i*zs-(zs-1):i*zs,find(auROCp==max(auROCp)));
    end
    phi_before_avg = mean(phi_before1z,1);

    phi_after1z = [];
    for i = 1:(size(h_after,1)/zs)
        phi_after1z(i,:) = h_after(i*zs-(zs-1):i*zs,find(auROCp==max(auROCp)));
    end
    phi_after_avg = mean(phi_after1z,1);
    subplot(3,1,3); plot([phi_before_avg';phi_after_avg']); title('average activity before & after lever press of highest auROC neuron');
    hold on; plot([zs+0.5,zs+0.5],[-0.3,0.5], '--');

    plo_before1z = [];
    for i = 1:(size(p_before,1)/zs)
        plo_before1z(i,:) = p_before(i*zs-(zs-1):i*zs,find(auROCp==min(auROCp)));
    end
    plo_before_avg = mean(plo_before1z,1);

    plo_after1z = [];
    for i = 1:(size(h_after,1)/zs)
        plo_after1z(i,:) = h_after(i*zs-(zs-1):i*zs,find(auROCp==min(auROCp)));
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
    for i = 1:size(z.dip_present,1)
        plot([z.dip_present(i,1),z.dip_present(i,1)],[-2,7], '--k', 'linewidth', 0.8);
    end

     subplot(3,1,2); plot(smoothed_rawdata(auROCd==min(auROCd),:)'); 
    title('lowest dipper presentation auROC neuron'); hold on; 
    for i = 1:size(z.dip_present,1)
        plot([z.dip_present(i,1),z.dip_present(i,1)],[-2,7], '--k', 'linewidth', 0.8);
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
    for i = 1:size(z.rewarded_entries,1)
        plot([z.rewarded_entries(i,1),z.rewarded_entries(i,1)],[-2,7], '--k', 'linewidth', 0.8);
    end

    subplot(3,1,2); plot(smoothed_rawdata(auROCr==min(auROCr),:)'); 
    title('lowest rewarded head entry auROC neuron'); hold on; 
    for i = 1:size(z.rewarded_entries,1)
        plot([z.rewarded_entries(i,1),z.rewarded_entries(i,1)],[-2,7], '--k', 'linewidth', 0.8);
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
    for i = 1:size(z.unrewarded_entries,1)
        plot([z.unrewarded_entries(i,1),z.unrewarded_entries(i,1)],[-2,7], '--k', 'linewidth', 0.8);
    end

    subplot(3,1,2); plot(smoothed_rawdata(auROCu==min(auROCu),:)'); 
    title('lowest unrewarded head entry auROC neuron'); hold on; 
    for i = 1:size(z.unrewarded_entries,1)
        plot([z.unrewarded_entries(i,1),z.unrewarded_entries(i,1)],[-2,7], '--k', 'linewidth', 0.8);
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
showclusters = @(neuron_image,indexhigh,indexlow) 2*max(neuron_image(:,:,:)>0,[],3) + max(neuron_image(:,:,indexhigh)>0,[],3)-max(neuron_image(:,:,indexlow)>0,[],3);

% Calculate top and bottom groups for all behaviors
% lever presses
top10_lp = top10(indexp); bottom10_lp = bottom10(indexp);
top10_r = top10(indexr); bottom10_r = bottom10(indexr);
top10_d = top10(indexd); bottom10_d = bottom10(indexd);
top10_u = top10(indexu); bottom10_u = bottom10(indexu);

% lever presses
if allfigures; figure; 
subplot(2,2,1); imagesc(showclusters(neuron_image,top10_lp,bottom10_lp)); title(sprintf('Top 10 for day %1.0f and time = %1.0f s\n \n Lever press neurons',day,zs/5));
% dipper presentations
subplot(2,2,2); imagesc(showclusters(neuron_image,top10_d,bottom10_d)); title(sprintf('[BLUE=inhibited; ORANGE=excited]\n \n Dipper presentation neurons'));
% rewarded entries
subplot(2,2,3); imagesc(showclusters(neuron_image,top10_r,bottom10_r)); title('Rewarded entry neurons');
% unrewarded entries
subplot(2,2,4); imagesc(showclusters(neuron_image,top10_u,bottom10_u)); title('Unrewarded entry neurons');

for i = 1:4; subplot(2,2,i); colormap(jet(5)); end
end
%%
% Outputs to SAVE to a master structure (called data_by_day)

data_by_day7677(day,zs).auROCp = auROCp;
data_by_day7677(day,zs).auROCd = auROCd;
data_by_day7677(day,zs).auROCu = auROCu;
data_by_day7677(day,zs).auROCr = auROCr;

data_by_day7677(day,zs).dr_5p=dr_5p;
data_by_day7677(day,zs).dr_95p=dr_95p;
data_by_day7677(day,zs).pr_5p=pr_5p;
data_by_day7677(day,zs).pr_95p=pr_95p;
data_by_day7677(day,zs).hr_5p=hr_5p;
data_by_day7677(day,zs).hr_95p=hr_95p;
data_by_day7677(day,zs).rr_5p=rr_5p;
data_by_day7677(day,zs).rr_95p=rr_95p;
data_by_day7677(day,zs).ur_5p=ur_5p;
data_by_day7677(day,zs).ur_95p=ur_95p;

data_by_day7677(day,zs).auROCpr_inc = auROCpr_inc;
data_by_day7677(day,zs).auROCdr_inc = auROCdr_inc;
data_by_day7677(day,zs).auROChr_inc = auROChr_inc;
data_by_day7677(day,zs).auROCrr_inc = auROCrr_inc;
data_by_day7677(day,zs).auROCur_inc = auROCur_inc;

data_by_day7677(day,zs).ci_outsideu = ci_outsideu;
data_by_day7677(day,zs).net_distu = net_distu;
data_by_day7677(day,zs).ci_outsider = ci_outsider;
data_by_day7677(day,zs).net_distr = net_distr;

data_by_day7677(day,zs).cell_rankr = cell_rankr;
data_by_day7677(day,zs).mean_indexr = mean_indexr; 
data_by_day7677(day,zs).dist_from_centerr = dist_from_centerr;
data_by_day7677(day,zs).rank_stdr = rank_stdr;
data_by_day7677(day,zs).centerr = centerr;
data_by_day7677(day,zs).correlationr = correlationr;
data_by_day7677(day,zs).cell_ranku = cell_ranku;
data_by_day7677(day,zs).mean_indexu = mean_indexu; 
data_by_day7677(day,zs).dist_from_centeru = dist_from_centeru;
data_by_day7677(day,zs).rank_stdu = rank_stdu;
data_by_day7677(day,zs).centeru = centeru;
data_by_day7677(day,zs).correlationu = correlationu;
data_by_day7677(day,zs).rank_imager = rank_imager;
data_by_day7677(day,zs).rank_imageu = rank_imageu;

data_by_day7677(day,zs).cluster_lp = showclusters(neuron_image,top10_lp,bottom10_lp);
data_by_day7677(day,zs).cluster_d = showclusters(neuron_image,top10_d,bottom10_d);
data_by_day7677(day,zs).cluster_r = showclusters(neuron_image,top10_r,bottom10_r);
data_by_day7677(day,zs).cluster_u = showclusters(neuron_image,top10_u,bottom10_u);

data_by_day7677(day,zs).corrcoef = corrcoef([auROCp(indexd);auROCd(indexd);auROCr(indexd);auROCu(indexd)]');

fprintf('Done with day %1.0f, zone size %1.0f \n',day,zs);