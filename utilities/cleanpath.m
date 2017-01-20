function [ output_path ] = cleanpath( input_path )
% cleanpath checks whether the entered path contains right separator and
% adds a separator at the end if it didn't have one
pos_filesep = {'/','\'};

% check separator
path_split = strsplit(input_path, pos_filesep);
output_path = [];
for ii = 1:length(path_split)
    if isempty(path_split{ii}) && ii == length(path_split), break; end
    output_path = [output_path path_split{ii} filesep];
end
if output_path(end) ~= filesep
    output_path = [output_path filesep];
end
end

