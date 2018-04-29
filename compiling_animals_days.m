%% Summarizes data from an rr5data structure
%
%   The rr5data structure should have:
%       -columns corresponding to mice
%       -rows corresponding to the days 
%
%   After the script is run, the structure rr5summary contains query-able
%   fields using 
%
%   rr5summary(# of day).mousename.one of the following:
%
%               head_entries 
%               dip_present 
%               lev_press 
%               inter_press_interval
%               rewarded_entries 
%               unrewarded_entries 
%               time_to_dipper 
%               dipper_received 
%
%  Any of these fields can be queried in the second cell of this code
%  simply set query = '' to any of the above (on line 89)

clear all
load('rr5data.mat')

days = size(rr5data,2);
animals = fields(rr5data);

for i = 1:size(animals,1)
    for j=1:days
        
        if ~evalin('base',sprintf('isempty(rr5data(%1.0f).%s)',j,animals{i}))
            clearvars tmp*

            assignin('base','tmp',evalin('base',sprintf('rr5data(%1.0f).%s',j,animals{i})))
            tmp(find(tmp(:,1)==0),:)=[]; % Eliminate random zeroes

            tmp_head = tmp(:,1)*5; %converts seconds to frames and caps at 1500 frames
            tmp_dipper = tmp(:,2)*5;
            tmp_press = tmp(:,3)*5;

            tmp_head_entries = tmp_head(find(tmp_head<1500));
            tmp_dip_present = tmp_dipper(find(tmp_dipper<1500));
            tmp_rlev_press = tmp_press(find(tmp_press<1500));

            % Getting the rewarded versus unrewarded entries
                rcounter = 1;
                ucounter = 1;
                for k = 1:size(tmp_head_entries,1)
                    thept = intersect(find(tmp_head_entries(k) > tmp_dip_present),find(tmp_head_entries(k) < (tmp_dip_present+25)));
                    if ~isempty(thept) % If there is a point in the range of 25 frames from the head entry where the dipper is up
                        tmp_rewarded_entries(rcounter,1) = tmp_head_entries(k);
                        tmp_dipper_received(rcounter,1) = tmp_dip_present(thept);
                        rcounter = rcounter + 1;
                    else
                        tmp_unrewarded_entries(ucounter,1) = tmp_head_entries(k);
                        ucounter = ucounter + 1;
                    end
                end

            % Time between dipper and head entry
            tmp_time_to_dipper = 0.2*(tmp_rewarded_entries-tmp_dipper_received);
            tmp_ipi = diff(tmp_rlev_press);

            evalin('base',sprintf('rr5summary(%1.0f).%s.%s=%s;',j,animals{i},'head_entries','tmp_head_entries'));
            evalin('base',sprintf('rr5summary(%1.0f).%s.%s=%s;',j,animals{i},'dip_present','tmp_dip_present'));
            evalin('base',sprintf('rr5summary(%1.0f).%s.%s=%s;',j,animals{i},'lev_press','tmp_rlev_press'));
            evalin('base',sprintf('rr5summary(%1.0f).%s.%s=%s;',j,animals{i},'inter_press_interval','tmp_ipi'));
            evalin('base',sprintf('rr5summary(%1.0f).%s.%s=%s;',j,animals{i},'rewarded_entries','tmp_rewarded_entries'));
            evalin('base',sprintf('rr5summary(%1.0f).%s.%s=%s;',j,animals{i},'unrewarded_entries','tmp_unrewarded_entries'));
            evalin('base',sprintf('rr5summary(%1.0f).%s.%s=%s;',j,animals{i},'time_to_dipper','tmp_time_to_dipper'));
            evalin('base',sprintf('rr5summary(%1.0f).%s.%s=%s;',j,animals{i},'dipper_received','tmp_dipper_received'));

            if numel(tmp_dip_present)~=numel(tmp_rewarded_entries)
                fprintf('A dipper was missed on day %1.0f by %s\n',j,animals{i});
            end
        end % End of if statement checking that there is data for that day
    end
end

clearvars -except rr5data rr5summary animals days


%% Plot all of the times to dipper
query_values_cell = {};
all_query_values = [];

query = 'inter_press_interval';

for i = 1:size(animals,1)
    for j=1:days

        if ~evalin('base',sprintf('isempty(rr5summary(%1.0f).%s)',j,animals{i}))
            
            assignin('base','tmp',evalin('base',sprintf('rr5summary(%1.0f).%s.%s',j,animals{i},query)));
            query_values_cell{i,j} = tmp; 
            all_query_values = [all_query_values; query_values_cell{i,j}];
            
        end
        
    end
end

figure; hist(all_query_values,300)
