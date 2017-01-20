function [sorter_table,sorter_names,nsorter,sorter_Varname,total_cond] = sorter_trials(data,sorterIn)
%sorter_trials generates sorter_table by selecting trials in the data table that pass each sorter condition
%   output:
%   sorter_table: array of logicals indicating trials in each condition [ntrial x ncond]
%   sorter_names: name of each sorter (e.g. CuingType)
%   nsorter: number of sorters
%   sorter_Varname: condition names (e.g. CuingType1, CuingType2)
%   total_cond: number of conditions

%   inputs:
%   data: data table ntrials x variables
%   sorterIn: structures containing sorter variables & values
%       ex.
%       CuingType: [1 2]
%       BiasType: [0 6]

trialVar_names = data.Properties.VariableNames;
ntrial = height(data);

sorter_names = fieldnames(sorterIn);
sorter_names = sorter_names'; % column to row vector
nsorter = length(sorter_names);
sorter_val = cell(size(sorter_names));
sorter_all = cell(size(sorter_names));
sorter_ind = cell(size(sorter_names));
ncond = zeros(size(sorter_names));
for st = 1:nsorter
    if any(ismember(trialVar_names,sorter_names{st}))
        sorter_val{st} = sorterIn.(sorter_names{st});
        ncond(st) = length(sorter_val{st});
        sorter_all{st} = zeros(ntrial,ncond(st));
        for cc = 1:ncond(st)
            sorter_all{st}(:,cc) = (data.(sorter_names{st}) == sorter_val{st}(cc));
        end
        sorter_ind{st} = 1:ncond(st);
    else
        disp([sorter_names{st} 'is not a valid sorter. No matched trial variables found']);
    end
end

% generate condition array (all combinations of sorters) then sorter table
tmp_cond = cell(size(sorter_names));
[tmp_cond{:}] = ndgrid(sorter_ind{:}); % sorter_val has to be a row vector e.g. {[1 2], [1 2], [4 5]}
cond_array = cell2mat(cellfun(@(v)v(:), tmp_cond, 'UniformOutput',false));
total_cond = size(cond_array,1);
sorter_table = ones(ntrial,total_cond); % sorter table starts from 1 then discard trials that don't pass each filter
sorter_Varname = cell(1,total_cond);
for cc = 1:total_cond
    for st = 1:nsorter % size(cond_array,2) must be equal to length(sorter_names)
        sorter_table(:,cc) = sorter_table(:,cc) & sorter_all{st}(:,cond_array(cc,st));
        
        if isempty(sorter_Varname{cc})
            sorter_Varname{cc} = [sorter_names{st} num2str(sorter_val{st}(cond_array(cc,st)))];
        else
            sorter_Varname{cc} = [sorter_Varname{cc} '_' sorter_names{st} num2str(sorter_val{st}(cond_array(cc,st)))];
        end
    end
end


end

