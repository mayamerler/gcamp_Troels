%% Reanalzying eduardos data with auROC
set(0,'defaultfigurecolor','w');
set(0,'defaultlinelinewidth',1)

rewards = m255.beh_time(find(m255.rewarded_beh==1))
norewards = m255.beh_time(find(m255.rewarded_beh==0))
lp = m255.press_time;

figure();

for i = 1:numel(rewards)
    t = rewards(i);
    aurocs_rewards(i,:) = mayaauroc(m255.Craw,[t-15:t],[t+1:t+15]);
    subplot(2,3,1); hold on; plot(sort(aurocs_rewards(i,:)),'o-','markerfacecolor','w','markersize',4); ylim([0,1]);
    xlim([0,105]);
    title('Rewarded','color','g');
end

for i = 1:numel(norewards)
    t = norewards(i);
    aurocs_norewards(i,:) = mayaauroc(m255.Craw,[t-15:t],[t+1:t+15]);
    subplot(2,3,2); hold on; plot(sort(aurocs_norewards(i,:)),'o-','markerfacecolor','w','markersize',4); ylim([0,1]);
    xlim([0,105]);
    title('Unrewarded');
end

for i = 1:numel(lp)
    t = lp(i);
    aurocs_lp(i,:) = mayaauroc(m255.Craw,[t-15:t],[t+1:t+15]);
    subplot(2,3,3); hold on; plot(sort(aurocs_lp(i,:)),'o-','markerfacecolor','w','markersize',4); ylim([0,1]);
    xlim([0,105]);
    title('Lever presses','color','r');
end

subplot(2,3,[4:6]); plot(sort(mean(aurocs_rewards,1)),'g'); hold on;
subplot(2,3,[4:6]); plot(sort(mean(aurocs_norewards,1)),'k');
subplot(2,3,[4:6]); plot(sort(mean(aurocs_lp,1)),'r');