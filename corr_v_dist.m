%% Loads data for mouse 8268 over 5 days

files{1} = 'D:\troels gcamp\Workspaces\Workspace 8268 RR5 Day1 calculated.mat';
files{2} = 'D:\troels gcamp\Workspaces\Workspace 8268 RR5 Day2 calculated.mat';
files{3} = 'D:\troels gcamp\Workspaces\Workspace 8268 RR5 Day3 calculated.mat';
files{4} = 'D:\troels gcamp\Workspaces\Workspace 8268 RR5 Day4 calculated.mat';
files{5} = 'D:\troels gcamp\Workspaces\Workspace 8268 RR5 Day5 calculated.mat';

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
day=2; % Set day to a value between 1 and 5, run the whole script to generate figures

load(files{day});
%load('D:\troels gcamp\Workspaces\rr5data.mat'); % This is to get the RR5 data without Excel 

zs=15; % sets size of time zone in samples (sampling rate for this data is 5 Hz

widegraph = [45 70 1000 550];

set(0,'defaultlinelinewidth',1.5);
set(0,'defaultfigurecolor','w');

%% find location of each neuron
neuron_image=[];

for i = 1:size(neuron.A,2)
    neuron_image(:,:,i) = reshape(neuron.A(:,i),[size(neuron.A,1)^(1/2),size(neuron.A,1)^(1/2)]);
end

n_image = double(sum(neuron_image,3)>0);


%% find correlation coefficients and distances between each set of neurons
clearvars neuron_CCs x_dist y_dist neuron_dist
neuron_CCs = corrcoef(neuron.C_raw');

for i = 1:size(neuron.C_raw,1)
    for j = 1:size(neuron.C_raw,1)
        x_dist = abs(center(i,1)-center(j,1));
        y_dist = abs(center(i,2)-center(j,2));
        neuron_dist(i,j) = (x_dist^2 + y_dist^2)^(1/2);
    end
end

m8268(day).neuron_CCs = neuron_CCs;
m8268(day).neuron_dist = neuron_dist;

%%

figure;

for i = 1:5
    plot(m7677(i).neuron_dist, m7677(i).neuron_CCs,'ob'); hold on;
    plot(m7675(i).neuron_dist, m7675(i).neuron_CCs,'or'); hold on;
end
for i = [1,3:5]
    plot(m8170(i).neuron_dist, m8170(i).neuron_CCs,'og'); hold on;
end
for i = 1:2
    plot(m8268(i).neuron_dist, m8268(i).neuron_CCs,'om'); hold on;
end

%% pool all days for each mouse
for i=1:5
    m7677_CCs = m7677(i).neuron_CCs((triu(m7677(i).neuron_CCs,1)~=0));
    m7675_CCs = m7675(i).neuron_CCs((triu(m7675(i).neuron_CCs,1)~=0));
    m7677_dist = m7677(i).neuron_dist((triu(m7677(i).neuron_dist,1)~=0));
    m7675_dist = m7675(i).neuron_dist((triu(m7675(i).neuron_dist,1)~=0));
end
for i=[1,3:5]
    m8170_CCs = m8170(i).neuron_CCs((triu(m8170(i).neuron_CCs,1)~=0));
    m8170_dist = m8170(i).neuron_dist((triu(m8170(i).neuron_dist,1)~=0));
end
for i = 1:2
    m8268_CCs = m8268(i).neuron_CCs((triu(m8268(i).neuron_CCs,1)~=0));
    m8268_dist = m8268(i).neuron_dist((triu(m8268(i).neuron_dist,1)~=0));
end


%% separate into distance categories
%all
[num_all,edges_all] = histcounts(m7677_CCs);
figure; plot([-0.475:0.05:0.55],num_all);
figure;
% separate by 25 pixel bins
for i = 1:10
    dist_CCs{i} = m7677_CCs(find(m7677_dist >= i*25-25 & m7677_dist <= i*25));
    [num{i},edges{i}] = histcounts(dist_CCs{i});
    xes{i} = [edges{i}(1) + ((edges{i}(2)-edges{i}(1))/2):((edges{i}(2)-edges{i}(1))/2)*2:edges{i}(size(edges{i},2))];
    plot(xes{i},num{i}); hold on;
end