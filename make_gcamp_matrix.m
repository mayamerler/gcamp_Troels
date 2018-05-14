% Parameters for analysis
% Set day to a value between 1 and 5, run the whole script to generate figures
mouse={'7674','7675','7677','8170','8265','8268'};

mydir = 'D:/troels gcamp/Workspaces';

for i = 1:numel(mouse)
    files = {};
    files{1} = sprintf('%s/Workspace %s RR5 Day1 calculated.mat',mydir,mouse{i});
    files{2} = sprintf('%s/Workspace %s RR5 Day2 calculated.mat',mydir,mouse{i});
    if strcmp(mouse{i},'8268')
    else
        files{3} = sprintf('%s/Workspace %s RR5 Day3 calculated.mat',mydir,mouse{i});
        files{4} = sprintf('%s/Workspace %s RR5 Day4 calculated.mat',mydir,mouse{i});
        if strcmp(mouse{i},'7677');
            files{5} = sprintf('%s/Workspace %s RR5 Day5 original.mat',mydir,mouse{i});
        else
            files{5} = sprintf('%s/Workspace %s RR5 Day5 calculated.mat',mydir,mouse{i});
        end
    end
    
for day = 1:size(files,2)
            
    load(files{day});
    tmp_neuron = neuron.C_raw;
    tmp_neu_image = zeros(size(neuron.A,1)^(1/2),size(neuron.A,1)^(1/2),size(neuron.A,2));
    
    for j = 1:size(neuron.A,2)
        tmp_neu_image(:,:,j) = reshape(neuron.A(:,j),[size(neuron.A,1)^(1/2),size(neuron.A,1)^(1/2)]);
    end
    
    evalin('base',sprintf('gcamp(%1.0f).m%s.%s=%s;',day,mouse{i},'neuron','tmp_neuron'));
    evalin('base',sprintf('gcamp(%1.0f).m%s.%s=%s;',day,mouse{i},'neu_image','tmp_neu_image'));
    
end

end