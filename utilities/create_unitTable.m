function [fin_table, overwrite] = create_unitTable(paths,recordingFile,overwrite)
%create_unitTable generates unit table sheets summarising recording units automatically
%   paths.recordingLog: path of recordingLog files
%   paths.recordingFile: path of xls files
%   paths.dataFile: path of plexon files
%   recordingFile: xls filename
%   overwrite: 0 - append unit_table to existing xls file | 1 - create new xls file

df_sheets = {'monkeys','dataFiles','chambers','units'};
df_header_monkeys = {'Monkey','ImplantID', 'Area', 'Stereotactic_AP', 'Stereotactic_ML', 'SurgeryDate', 'Comment'};
df_header_dataFiles = {'Date','Rec', 'DataFile', 'TypeOfFile', 'Monkey', 'Task', 'TaskVariation', 'Experimenter', 'Setup', 'Comment'};
df_header_chambers = {'Date', 'DataFile', 'Area', 'Drive','Head', 'Matrix_x','Matrix_y', 'GuideTubeDepth', 'GuideTubeLength', 'nElectrode', 'nChannel', 'nUnit', 'Comment'};
df_header_units = {'Date', 'Rec', 'DataFile', 'Monkey', 'Area', 'Task', 'TaskVariation', 'Channel', 'Unit', 'LFPch', 'LPFcomment','depth', 'Start','Stop', 'Signal', 'Isolation', 'Comment', 'SortMethod', 'Sorter','TypeOfFile'};

df_area = {'PMd';'PRR'};
df_electrodeInChambers = [1:5;6:10]; % default: minimatrix PMd:1-5 | PRR:6-10

df_nSheet = length(df_sheets);
df_nField = [length(df_header_monkeys),length(df_header_dataFiles),length(df_header_chambers),length(df_header_units)];

% some limits
% (this setting is used for 2 minimatrix recordings: 5x2 channels)
max_unit = 10; % max. number of units in 1 channel
max_chan = 10; % max. number of recording channels
chan_digit = numel(num2str(max_chan)); % used for unit name
possible_unit_id = 1:max_unit; % start at 1 - ignore unit 0 (unsorted)
unit_suffix = cellstr(char('a' + (1:max_unit)-1)')';

if ~exist(paths.recordingLog,'dir')
    paths.recordingLog = uigetdir(paths.dataFile, 'Select a folder containing recording log files');
end
paths.recordingLog = cleanpath(paths.recordingLog);

if ~exist(paths.recordingFile,'dir')
    paths.recordingFile = uigetdir('', 'Select a folder containing unit table');
end
paths.recordingFile = cleanpath(paths.recordingFile);


% check whether the file already exists
recordingFullFile = [paths.recordingFile recordingFile];
unitTable_exist = 0;
if exist(recordingFullFile,'file')
    unitTable_exist = 1;
    unitTable_format = recordingFullFile(end-3:end);
elseif exist([recordingFullFile '.xls'],'file')
    unitTable_exist = 1;
    unitTable_format = '.xls';
    recordingFullFile = [recordingFullFile '.xls'];
elseif exist([recordingFullFile '.csv'],'file')
    unitTable_exist = 1;
    unitTable_format = '.csv';
    recordingFullFile = [recordingFullFile '.csv'];
end

if unitTable_exist
    recordingFullFile_old = [paths.recordingFile 'old_' datestr(now,'yymmddHHMM') '_' recordingFile];
    if overwrite
        movefile(recordingFullFile,recordingFullFile_old);
        disp(['existing file "' recordingFile '" will be overwritten...'])
        disp(['information in the old file is moved to "' recordingFullFile_old '"'])
    else
        if strcmp(unitTable_format,'.xls')
            [~, old_sheets] = xlsfinfo(recordingFullFile);
            old_nSheet = length(old_sheets);
        elseif strcmp(unitTable_format,'.csv')
            old_sheets = [];  % TODO: I have no csv yet
            old_nSheet = 0;   % TODO
        end
        
        if ~all(ismember(df_sheets,old_sheets))  % if existing file doesn't contain all default sheets -> rename the old one and creat a new xls file
            movefile(recordingFullFile,recordingFullFile_old);
            overwrite = 1;
            disp(['existing file "' recordingFile '" does not contain all default tables...'])
            disp(['A new file is created. Information in the old file is moved to "' recordingFullFile_old '"'])
        else
            old_table = cell(1,old_nSheet);
            for sh = 1:old_nSheet
                old_table{sh} = readtable(recordingFullFile,'Sheet',df_sheets{sh});
                if width(old_table{sh}) ~= df_nField(sh)   % if existing file doesn't contain all default fields -> rename the old one and creat a new xls file
                    movefile(recordingFullFile,recordingFullFile_old);
                    overwrite = 1;
                    disp(['existing file "' recordingFile '" does not contain all default fields...'])
                    disp(['A new file is created. Information in the old file is moved to "' recordingFullFile_old '"'])
                    break;
                end
            end
        end
    end
else
    overwrite = 1;
end

listing_files = dir([paths.dataFile,'*.plx']);
dataFiles = cell(length(listing_files),1);
for ii = 1:length(listing_files)
    dataFiles{ii} = listing_files(ii).name;
end

nFile = length(dataFiles);

typeOfFile = cell(nFile,1);
for ii = 1:nFile
    typeOfFile{ii} = dataFiles{ii}(end-3:end);
end

% to be tested when we have header files for all recording files
% look for recordingLog files
useLogFiles = 0;
listing_log = dir([paths.recordingLog,'*recordingLog.mat']);
nLogFile = size(listing_log,1);
if nLogFile
    lognames_prefix = cell(nLogFile,1);
    for ii = 1:nLogFile
        lognames_prefix{ii} = listing_log(ii).name(1:10);
    end
end

% variables for 'dataFiles' table
date_2        = cell(nFile,1);
rec_2         = cell(nFile,1);
dataFile_2    = cell(nFile,1);
monkey_2      = cell(nFile,1);
experimenter  = cell(nFile,1);
setup         = cell(nFile,1);
matrix        = cell(nFile,1);                 
task_2        = cell(nFile,1);
task_variation_2 = cell(nFile,1);
comment_2 = cell(nFile,1);


% variables for 'chambers' table
date_3 = [];
dataFile_3 = [];
area_3 = [];
drive = [];
head = [];
matrix_x = [];
matrix_y = [];
ZDepth = [];
guideTube = [];
nElectrode = [];
nChannel = [];
nUnit = [];
comment_3 = [];

% variables for 'units' table
date_4 = [];
dataFile_4 = [];
rec_4 = [];
monkey_4 = [];
area_4 = [];
task_4 = [];
task_variation_4 = [];
channel = [];
unit = [];
depth = [];
LFPch = [];

plxUnits = cell(nFile,2);
plxUnitsSum = cell(nFile,1);
for ii = 1:nFile
    dataFile_2{ii} = dataFiles{ii}(1:end-4);
    data_prefix = dataFiles{ii}(1:10);
    useLogFiles = 0; 
    if nLogFile
        if any(ismember(lognames_prefix,data_prefix))
            useLogFiles = 1;
            ind_log = find(ismember(lognames_prefix,data_prefix),1,'last');
        end
    end
  
    % take the 1st level header info from recordingLog
    if useLogFiles
        load([paths.recordingLog listing_log(ind_log).name])
        date_2{ii} = datestr(header_var.date,'dd.mm.yyyy'); 
        monkey_2{ii} = header_var.monkey(1:3);
        task_2{ii} = header_var.task;
        tmp_nChamber = size(header_var.recordingAreas{1},1);
        switch tmp_nChamber
            case 1, matrix{ii} = 'SM';
            case 2, matrix{ii} = 'DM';
        end
        % number of recording has to come from the filename
        rec_2{ii} = sscanf(dataFiles{ii},['%*11s' '%2s' '*']);
        task_variation_2{ii} = [];
        experimenter{ii} = header_var.recorder;
        setup{ii} = header_var.setup;
    else
        % if not extract data from filename: given the filename obeys a centain format
        tmp_splitName = textscan(dataFiles{ii},'%s %s %s %s %s %s %s','Delimiter','_');
        jj = size(tmp_splitName,2);
        while jj >= 1
            if isempty(tmp_splitName{jj})
                tmp_splitName{1,jj} = [];
            elseif strfind(tmp_splitName{jj}{1},'.plx')
                tmp_splitName{1,jj} = [];
            else
                tmp_splitName{1,jj} = tmp_splitName{jj}{1};
            end
            jj = jj - 1;
        end
        monkey_2{ii} = sscanf(tmp_splitName{1},['%3s' '*']);                          % first 3 letters : monkey's name
        tmp_date = datenum(sscanf(tmp_splitName{1},['%*4s' '%6s']),'yyMMdd');       % next 6 letters : date
        date_2{ii} = datestr(tmp_date,'dd.MM.yyyy');                                  % convert to xls format
        rec_2{ii} = tmp_splitName{2};                                                 % next 2 letters : recording number
        matrix{ii} = tmp_splitName{3};                                              % next 2 letters : single or double matrix recording
        task_2{ii} = tmp_splitName{4};                                                % next letters : task name
        task_variation_2{ii} = tmp_splitName{5};                                      % next letters : task variations
        switch matrix{ii}
            case 'SM', tmp_nChamber = 1;
            case 'DM', tmp_nChamber = 2;
        end
        experimenter{ii} = [];
        setup{ii} = [];
    end
    
    date_3 = [date_3;repmat(date_2(ii),tmp_nChamber,1)];
    dataFile_3 = [dataFile_3;repmat(dataFile_2(ii),tmp_nChamber,1)];
    comment_3 = [comment_3;cell(tmp_nChamber,1)];
    if useLogFiles
        area_3 = [area_3;subHeader_var.Properties.RowNames];
        drive = [drive;subHeader_var.matrixNumber];
        head = [head;subHeader_var.matrixHead];
        matrix_x = [matrix_x;num2cell(subHeader_var.guideTubeX)];
        matrix_y = [matrix_y;num2cell(subHeader_var.guideTubeY)];
        % add a space in front of guideTubeZ to avoid autoconversion to date in xls
        for zz = 1:size(subHeader_var.guideTubeZ,1)
            subHeader_var.guideTubeZ{zz} = [' ' subHeader_var.guideTubeZ{zz}];
        end
        ZDepth = [ZDepth;subHeader_var.guideTubeZ];
        guideTube = [guideTube;num2cell(subHeader_var.guideTubeLength)];
        tmp_nE = zeros(tmp_nChamber,1);
        for cc = 1:tmp_nChamber
            tmp_nE(cc) = length(subHeader_var.channels{cc});
        end
        nElectrode = [nElectrode;num2cell(tmp_nE)];
    else
        drive = [drive;cell(tmp_nChamber,1)];
        head = [head;cell(tmp_nChamber,1)];
        matrix_x = [matrix_x;cell(tmp_nChamber,1)];
        matrix_y = [matrix_y;cell(tmp_nChamber,1)];
        ZDepth = [ZDepth;cell(tmp_nChamber,1)];
        guideTube = [guideTube;cell(tmp_nChamber,1)];
        nElectrode = [nElectrode;cell(tmp_nChamber,1)];
    end
    
    % information for each unit -> the 3rd sheet
    % count number of channels & units in each plx file
    
    filename = [paths.dataFile dataFiles{ii}];
    [plx.counts.tscounts, plx.counts.wfcounts, plx.counts.evcounts] = plx_info(filename,1);
    % tscounts, wfcounts are indexed by (unit+1,channel+1)
    % tscounts(:,ch+1) is the per-unit counts for channel ch
    % sum( tscounts(:,ch+1) ) is the total wfs for channel ch (all
    % units)
    
    [file_unit_id,file_channel_id] = find(plx.counts.tscounts);
    file_channel_id = file_channel_id-1;
    file_unit_id = file_unit_id-1;
    unit_id_to_extract = file_unit_id(file_unit_id>0); % discard unsorted unit (unit == 0)
    channel_id_to_extract = file_channel_id(file_unit_id>0);
    tmp_nU = size(channel_id_to_extract,1);
    
    date_4 = [date_4;repmat(date_2(ii),tmp_nU,1)];
    dataFile_4 = [dataFile_4;repmat(dataFile_2(ii),tmp_nU,1)];
    rec_4 = [rec_4;repmat(rec_2(ii),tmp_nU,1)];
    monkey_4 = [monkey_4;repmat(monkey_2(ii),tmp_nU,1)];
    task_4 = [task_4;repmat(task_2(ii),tmp_nU,1)];
    task_variation_4 = [task_variation_4;repmat(task_variation_2(ii),tmp_nU,1)];
    
    channel = [channel;num2cell(channel_id_to_extract)];
    
    tmp_nUnit = zeros(tmp_nChamber,1);
    tmp_nChannel = zeros(tmp_nChamber,1);
    tmp_area = [];
    
    if useLogFiles
        for cc = 1:tmp_nChamber
            if any(ismember(subHeader_var.channels{cc},channel_id_to_extract))
                tmp_nUnit(cc) = sum(ismember(channel_id_to_extract,subHeader_var.channels{cc}));
                tmp_nChannel(cc) = sum(ismember(subHeader_var.channels{cc},channel_id_to_extract));
                tmp_area = [tmp_area;repmat(subHeader_var.Properties.RowNames(cc),tmp_nUnit(cc),1)];
            end
        end
    else
        switch tmp_nChamber
            case 1
                for cc = 1:size(df_area,1)
                    if any(ismember(df_electrodeInChambers(cc,:),channel_id_to_extract))
                        tmp_nUnit = sum(ismember(channel_id_to_extract,df_electrodeInChambers(cc,:)));
                        tmp_nChannel = sum(ismember(df_electrodeInChambers(cc,:),channel_id_to_extract));
                        tmp_area = [tmp_area;repmat(df_area(cc),tmp_nUnit,1)];
                    end
                end
            case 2
                for cc = 1:tmp_nChamber
                    if any(ismember(df_electrodeInChambers(cc,:),channel_id_to_extract))
                        tmp_nUnit(cc) = sum(ismember(channel_id_to_extract,df_electrodeInChambers(cc,:)));
                        tmp_nChannel(cc) = sum(ismember(df_electrodeInChambers(cc,:),channel_id_to_extract));
                        tmp_area = [tmp_area;repmat(df_area(cc),tmp_nUnit(cc),1)];
                    end
                end
        end
    end
    area_4 = [area_4;tmp_area];
    
    for jj = 1:tmp_nU
        channel_id = channel_id_to_extract(jj);
        unit_id = unit_id_to_extract(jj);
        unit_name = [num2str(channel_id,['%0' num2str(chan_digit) 'u']) unit_suffix{possible_unit_id==unit_id}];
        unit = [unit;{unit_name}];
        tmp_lfp = {[]};
        LFPch = [LFPch;tmp_lfp];
    end
    
    % 3rd level data from single unit data
    nChannel = [nChannel;num2cell(tmp_nChannel)];
    nUnit = [nUnit;num2cell(tmp_nUnit)];
    if ~useLogFiles % check again
        switch tmp_nChamber
            case 1, area_3 = [area_3;unique(tmp_area)];
            case 2, area_3 = [area_3;df_area];
        end
    end
end

% variables for 'monkeys' table
monkey_list = unique(monkey_4);
nMonkey = length(monkey_list);
area_list = unique(area_4);
nArea = length(area_list);

monkey_1    = cell(nMonkey*nArea,1);
implantID   = cell(nMonkey*nArea,1);
area_1    = cell(nMonkey*nArea,1);
stereotactic_ap = cell(nMonkey*nArea,1);
stereotactic_ml = cell(nMonkey*nArea,1);
surgery_date = cell(nMonkey*nArea,1);                
comment_1 = cell(nMonkey*nArea,1);

for mm = 1:nMonkey
    for aa = 1:nArea
        monkey_1{((mm-1)*nArea)+aa} = monkey_list{mm};
        area_1{((mm-1)*nArea)+aa} = area_list{aa};
    end
end

empty_field = LFPch;
    
% write xls file (sheet1: sessions, sheet2: chambers, sheet3: units)
new_table = cell(1,df_nSheet);
fin_table = cell(1,df_nSheet);

new_table{1} = cell2table([monkey_1,implantID,area_1,stereotactic_ap,stereotactic_ml,surgery_date,comment_1],'VariableNames',df_header_monkeys);
new_table{2} = cell2table([date_2,rec_2,dataFile_2,typeOfFile,monkey_2,task_2,task_variation_2,experimenter,setup,comment_2],'VariableNames',df_header_dataFiles);
new_table{3} = cell2table([date_3,dataFile_3,area_3,drive,head,matrix_x,matrix_y,ZDepth,guideTube,nElectrode,nChannel,nUnit,comment_3],'VariableNames',df_header_chambers);
new_table{4} = cell2table([date_4,rec_4,dataFile_4,monkey_4,area_4,task_4,task_variation_4,channel,unit,LFPch,empty_field,empty_field,empty_field,empty_field,empty_field,empty_field,empty_field,empty_field,empty_field,empty_field],'VariableNames',df_header_units);

if ~overwrite
    % append new data files to existing table
    
    for sh = 1:df_nSheet
        % set recording numbers to be string
        if ismember('Rec',old_table{sh}.Properties.VariableNames), old_table{sh}.Rec = cellstr(num2str(old_table{sh}.Rec,'%02d')); end

        % check missing values and set to be empty cells
        MissingCol = old_table{sh}.Properties.VariableNames(all(ismissing(old_table{sh}),1));
        nMissingCol = length(MissingCol);
        if nMissingCol
            for nm = 1:nMissingCol
                old_table{sh}.(MissingCol{nm}) = num2cell(old_table{sh}.(MissingCol{nm}));
            end
        end
        % todo:check string and cellstr them
    end
    
    % monkey sheet
    new_monkey = setdiff(new_table{1}.Monkey,old_table{1}.Monkey);
    new_table{1}.Monkey = categorical(new_table{1}.Monkey);
    fin_table{1} = old_table{1};
    for nn = 1:length(new_monkey)
        fin_table{1} = [fin_table{1};new_table{1}(new_table{1}.Monkey == new_monkey(nn),:)];
    end
    fin_table{1} = sortrows(fin_table{1},{'Monkey','Area'});
    
    % all other sheets
    new_datafile = setdiff(new_table{2}.DataFile,old_table{2}.DataFile);
    for sh = 2:df_nSheet
        new_table{sh}.DataFile = categorical(new_table{sh}.DataFile);
        fin_table{sh} = old_table{sh};
        for nn = 1:length(new_datafile)
            fin_table{sh} = [fin_table{sh};new_table{sh}(new_table{sh}.DataFile == new_datafile(nn),:)];
        end
    end
    fin_table{4} = sortrows(fin_table{4},{'DataFile','Unit'});
    fin_table{3} = sortrows(fin_table{3},{'DataFile','Area'});
    fin_table{2} = sortrows(fin_table{2},'DataFile');
else
    for sh = 1:df_nSheet
        fin_table{sh} = new_table{sh};
    end
end