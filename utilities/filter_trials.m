function [filter_table,filter_names,nfilter] = filter_trials(data,filterIn)
%filter_trials generates filter_table by selecting trials in the data table that pass filter options
%   output:
%   filter_table: array of logicals indicating trials that pass the filter

%   inputs:
%   data: data table ntrials x variables
%   filterIn: structures containing filter variables & values
%       ex.
%       Success: 1
%       CuingType: 1

    trialVar_names = data.Properties.VariableNames;
    ntrial = height(data);
    
    filter_names = fieldnames(filterIn);
    nfilter = length(filter_names);
    filter_val = cell(size(filter_names));
    filter_table = ones(ntrial,1); % filter array starts from 1 then discard trials that don't pass each filter
    for ft = 1:nfilter
        if any(ismember(trialVar_names,filter_names{ft}))
            filter_val{ft} = filterIn.(filter_names{ft});
            tmp_filter = zeros(ntrial,1);
            for fv = 1:length(filter_val{ft})
                tmp_filter = tmp_filter | (data.(filter_names{ft}) == filter_val{ft}(fv));
            end
            filter_table = filter_table & tmp_filter;
        else
            disp([filter_names{ft} 'is not a valid filter. No matched trial variables found']);
        end
    end


end

