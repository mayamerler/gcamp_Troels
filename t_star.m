%% find auROC at each time point for each neuron
zs = 15;

for j = 1:size(neuron.C_raw,1)
    for i = 1:(1500-2*zs)
        zone1 = [i:i+zs];
        zone2 = [i+zs+1:i+2*zs];
        neu = neuron.C_raw(j,:);
        neuronROC(i+zs,j) = mayaauroc(neu, zone1, zone2);
    end
end

%% plot auROC over time with rewarded head entries

figure; plot(neuronROC(:,1), 'k'); title('auROC over time with rewarded head entries(r)'); hold on;
for h = 1:size(z.rewarded_entries,1)
    plot([z.rewarded_entries(h,1),z.rewarded_entries(h,1)],[0,1.0], '--r','linewidth',0.8);
end

%% distribution of time between dipper presentation and rewarded head entry

for i = 1:size(z.rewarded_entries,1)
    hd_int(i) = z.rewarded_entries(i)-z.dip_present(i);
end

%% convert all auROC points below 0.5 to the same distance above 0.5 (all troughs become peaks)

for i = 1:size(neuronROC,2)
    neu_abs(:,i) = abs(0.5-neuronROC(:,i))+0.5;
end
%% find local maxima in auROC after each rewarded head entry

for j = 1:size(neuronROC,2)
    for i = 1:size(z.rewarded_entries,1)
        [peakheightr(i,j),peaklocr(i,j)] = findpeaks(neu_abs((z.rewarded_entries(i)-1):z.rewarded_entries(i)+20,j), 'NPeaks', 1);
    end
end

figure; subplot(2,4,1); histogram(peaklocr); title('all t* rewarded entries');
subplot(2,4,5); histogram(peaklocr(peakheightr>0.8)); title('t* for peaks greater than 0.8');



% find local maxima in auROC after each dipper presentation

for j = 1:size(neuronROC,2)
    for i = 1:size(z.dip_present,1)
        [peakheightd(i,j),peaklocd(i,j)] = findpeaks(neu_abs(z.dip_present(i):end,j), 'NPeaks', 1);
    end
end

subplot(2,4,2); histogram(peaklocd); title('all t* dipper presentations');
subplot(2,4,6); histogram(peaklocd(peakheightd>0.8)); title('t* for peaks greater than 0.8');

% find local maxima in auROC after each unrewarded head entry

for j = 1:size(neu_abs,2)
    for i = 1:size(z.unrewarded_entries,1)-1
        
        [peakheightu(i,j),peaklocu(i,j)] = findpeaks(neu_abs(z.unrewarded_entries(i):end,j), 'NPeaks', 1);
    end
end

subplot(2,4,3); histogram(peaklocu); title('all t* unrewarded entries');
subplot(2,4,7); histogram(peaklocu(peakheightu>0.8)); title('t* for peaks greater than 0.8');

% find local maxima in auROC after each lever press

for j = 1:size(neuronROC,2)
    for i = 1:size(z.rlev_press,1)
        [peakheightp(i,j),peaklocp(i,j)] = findpeaks(neu_abs(z.rlev_press(i):end,j), 'NPeaks', 1);
    end
end

subplot(2,4,4); histogram(peaklocp); title('all t* lever presses');
subplot(2,4,8); histogram(peaklocp(peakheightp>0.8)); title('t* for peaks greater than 0.8');set(gcf,'position',widebox);
name = sprintf('D:/troels gcamp/cell rank figs/7675day%1.0f.svg',day);
saveas(gcf,name);

figure; plot(peaklocr,peakheightr,'o');




%% find local maxima in auROC after before rewarded head entry

for j = 1:size(neuronROC,2)
    for i = 1:size(z.rewarded_entries,1)
        [peakheightrb(i,j),peaklocrb(i,j)] = findpeaks(neu_abs((fix(z.rewarded_entries(i)+1)):-1:1,j), 'NPeaks', 1);
    end
end

figure; subplot(2,4,1); histogram(peaklocrb); title('all t* rewarded entries');
subplot(2,4,5); histogram(peaklocrb(peakheightrb>0.8)); title('t* for peaks greater than 0.8');

%%
for j = 1:size(neuronROC,2)
    for i = 1:size(z.rewarded_entries,1)
        [peakheightrb(i,j),peaklocrb(i,j)] = findpeaks(neu_abs(fix(z.rewarded_entries(i)-30):fix(z.rewarded_entries(i)),j), 'NPeaks', 1);
    end
end

figure; subplot(2,4,1); histogram(peaklocrb); title('all t* rewarded entries');
subplot(2,4,5); histogram(peaklocrb(peakheightrb>0.8)); title('t* for peaks greater than 0.8');



%%
% find local maxima in auROC after each dipper presentation

for j = 1:size(neuronROC,2)
    for i = 1:size(z.dip_present,1)
        [peakheightdb(i,j),peaklocdb(i,j)] = findpeaks(neu_abs(fix(z.dip_present(i)+1):-1:1,j), 'NPeaks', 1);
    end
end

subplot(2,4,2); histogram(peaklocdb); title('all t* dipper presentations');
subplot(2,4,6); histogram(peaklocdb(peakheightdb>0.8)); title('t* for peaks greater than 0.8');

% find local maxima in auROC after each unrewarded head entry

for j = 1:size(neu_abs,2)
    for i = 1:size(z.unrewarded_entries,1)-1
        
        [peakheightub(i,j),peaklocub(i,j)] = findpeaks(neu_abs((fix(z.unrewarded_entries(i)+1)):-1:1,j), 'NPeaks', 1);
    end
end

subplot(2,4,3); histogram(peaklocub); title('all t* unrewarded entries');
subplot(2,4,7); histogram(peaklocub(peakheightub>0.8)); title('t* for peaks greater than 0.8');

% find local maxima in auROC after each lever press

for j = 1:size(neuronROC,2)
    for i = 1:size(z.rlev_press,1)
        [peakheightpb(i,j),peaklocpb(i,j)] = findpeaks(neu_abs(fix(z.rlev_press(i)+1):-1:1,j), 'NPeaks', 1);
    end
end

subplot(2,4,4); histogram(peaklocpb); title('all t* lever presses');
subplot(2,4,8); histogram(peaklocpb(peakheightpb>0.8)); title('t* for peaks greater than 0.8');set(gcf,'position',widebox);
name = sprintf('D:/troels gcamp/cell rank figs/7675day%1.0f.svg',day);
saveas(gcf,name);