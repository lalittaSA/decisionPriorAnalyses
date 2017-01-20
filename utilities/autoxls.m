function autoxls(path,filenames,options)
%autoxls generates xls sheet summarising recording automatically
%   options.xxx
%   options.yyy

df_header_dataFiles = {'Date','Rec', 'DataFile', 'TypeOfFile', 'Monkey', 'Task', 'TaskVariation', 'GeneralComments'};
df_header_chambers = {'Date', 'DataFile', 'Area', 'Drive','Head', 'MatrixCoordinates', 'GuideTubeDepth', 'GuideTubeLength', 'nElectrode', 'nChannel', 'nUnit', 'Comment'};
df_header_units = {'Date', 'Rec', 'DataFile', 'Monkey', 'Area', 'Task', 'TaskVariation', 'Channel', 'Unit', 'LFPch', 'LPFcomment','depth', 'timeWindow', 'Signal', 'Isolation', 'Comment', 'SortMethod', 'Sorter','TypeOfFile'};

df_electrodeInChambers = [1:5;6:10]; % default: PMd:1-5 | PRR:6-10
df_area = {'PMd';'PRR'};

nSheet = 3;
% check whether the file already exists
recordingFullFile = [options.recordingPath options.recordingFile];
if exist(recordingFullFile,'file')
    if ~options.overwrite
%         [numData, txtData, rawData] = xlsread(recordingFullFile, 'units');
        old_table = cell(1,nSheet);
        old_table{1} = readtable(recordingFullFile,'Sheet','dataFiles');
        old_table{2} = readtable(recordingFullFile,'Sheet','chambers');
        old_table{3} = readtable(recordingFullFile,'Sheet','units');
%         xlsHeaderOS = xlsTxtOS(1,:);

        for sh = 1:nSheet
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
            
            % check string and cellstr them
            
        end
    else
        disp(['existing xls file "' options.recordingFile '" will be overwritten...'])
    end
else
    options.overwrite = 1;
end

nFile = length(filenames);

% get dataPaths
dataPath = repmat({path},nFile,1);

typeOfFile = cell(nFile,1);
for ii = 1:nFile
    typeOfFile{ii} = filenames{ii}(end-3:end);
end

% to be tested when we have header files for all recording files
% look for recordingLog files
useLogFiles = 0;
listing_log = dir([options.recordingLog,'\*.mat']);
nLogFile = size(listing_log,1);
if nLogFile
    lognames_prefix = cell(nLogFile,1);
    for ii = 1:nLogFile
        lognames_prefix{ii} = listing_log(ii).name(1:10);
    end
end

% variables for 'dataFiles' table
date        = cell(nFile,1);
dataFile    = cell(nFile,1);
rec         = cell(nFile,1);
monkey      = cell(nFile,1);
matrix      = cell(nFile,1);                 
task        = cell(nFile,1);
task_variation = cell(nFile,1);
general_comment = cell(nFile,1);

% variables for 'chambers' table
date_2 = [];
dataFile_2 = [];
area = [];
drive = [];
head = [];
% Xpos = [];
% Ypos = [];
matrixCoor = [];
ZDepth = [];
guideTube = [];
nElectrode = [];
nChannel = [];
nUnit = [];
comment = [];

% variables for 'units' table
date_3 = [];
dataFile_3 = [];
rec_3 = [];
monkey_3 = [];
area_3 = [];
task_3 = [];
task_variation_3 = [];
channel = [];
unit = [];
depth = [];
LFPch = [];

plxUnits = cell(nFile,2);
plxUnitsSum = cell(nFile,1);
for ii = 1:nFile
    dataNameSpk = [dataPath{ii},filenames{ii},'Spk','.mat'];
    dataNameAna = [dataPath{ii},filenames{ii},'Ana','.mat'];
    dataFile{ii} = filenames{ii}(1:end-4);
    data_prefix = filenames{ii}(1:10);
    useLogFiles = 0; 
    if nLogFile
        if any(ismember(lognames_prefix,data_prefix))
            useLogFiles = 1;
            ind_log = find(ismember(lognames_prefix,data_prefix),1,'last');
        end
    end
    
    % take the 1st level header info from recordingLog
    if useLogFiles
        load([options.recordingLog listing_log(ind_log).name])
        date{ii} = datestr(header_var.date,'dd.mm.yyyy'); 
        monkey{ii} = header_var.monkey(1:3);
        task{ii} = header_var.task;
        tmp_nChamber = size(header_var.recordingAreas{1},1);
        switch tmp_nChamber
            case 1, matrix{ii} = 'SM';
            case 2, matrix{ii} = 'DM';
        end
        % number of recording has to come from the filename
        rec{ii} = sscanf(filenames{ii},['%*11s' '%2s' '*']);
        task_variation{ii} = [];
    else
        % if not extract data from filename: given the filename obeys a centain format
        tmp_splitName = textscan(filenames{ii},'%s %s %s %s %s %s %s','Delimiter','_');
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
        monkey{ii} = sscanf(tmp_splitName{1},['%3s' '*']);                          % first 3 letters : monkey's name
        tmp_date = datenum(sscanf(tmp_splitName{1},['%*4s' '%6s']),'yyMMdd');       % next 6 letters : date
        date{ii} = datestr(tmp_date,'dd.MM.yyyy');                                  % convert to xls format
        rec{ii} = tmp_splitName{2};                                                 % next 2 letters : recording number
        matrix{ii} = tmp_splitName{3};                                              % next 2 letters : single or double matrix recording
        task{ii} = tmp_splitName{4};                                                % next letters : task name
        task_variation{ii} = tmp_splitName{5};                                      % next letters : task variations
        switch matrix{ii}
            case 'SM', tmp_nChamber = 1;
            case 'DM', tmp_nChamber = 2;
        end
    end
    % recorder? setup? task variation?
    
    date_2 = [date_2;repmat(date(ii),tmp_nChamber,1)];
    dataFile_2 = [dataFile_2;repmat(dataFile(ii),tmp_nChamber,1)];
    comment = [comment;cell(tmp_nChamber,1)];
    if useLogFiles
        area = [area;subHeader_var.Properties.RowNames];
        drive = [drive;subHeader_var.matrixNumber];
        head = [head;subHeader_var.matrixHead];
%         Xpos = [Xpos;num2cell(subHeader_var.guideTubeX)];
%         Ypos = [Ypos;num2cell(subHeader_var.guideTubeY)];
        tmp_coor = cell(tmp_nChamber,1);
        for xx = 1:size(subHeader_var.guideTubeX,1)
            tmp_coor{xx} = ['[' num2str(subHeader_var.guideTubeX(xx)) ',' num2str(subHeader_var.guideTubeY(xx)) ']'];
        end
        matrixCoor = [matrixCoor;tmp_coor];
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
%         Xpos = [Xpos;cell(tmp_nChamber,1)];
%         Ypos = [Ypos;cell(tmp_nChamber,1)];
        matrixCoor = [matrixCoor;cell(tmp_nChamber,1)];
        ZDepth = [ZDepth;cell(tmp_nChamber,1)];
        guideTube = [guideTube;cell(tmp_nChamber,1)];
        nElectrode = [nElectrode;cell(tmp_nChamber,1)];
    end
    
    % information for each unit -> the 3rd sheet
    % count number of channels & units in each plx file
    load(dataNameSpk)
    fn = fieldnames(data);
    nU = length(fn);
    
    date_3 = [date_3;repmat(date(ii),nU,1)];
    dataFile_3 = [dataFile_3;repmat(dataFile(ii),nU,1)];
    rec_3 = [rec_3;repmat(rec(ii),nU,1)];
    monkey_3 = [monkey_3;repmat(monkey(ii),nU,1)];
    task_3 = [task_3;repmat(task(ii),nU,1)];
    task_variation_3 = [task_variation_3;repmat(task_variation(ii),nU,1)];
    
    plxUnits{ii,2} = cell(nU,2);
    tmp_nUnit = zeros(tmp_nChamber,1);
    tmp_nChannel = zeros(tmp_nChamber,1);
    new_ch = 1;
    old_ch = 0;
    tmp_area = [];
    for jj = 1:nU
        % get channel numbers
        temp_ch0 = fn{jj}(end-1);
        temp_ch1 = fn{jj}(end-2);
        if ~isnan(str2double(temp_ch1))
            temp_ch = str2double([temp_ch1,temp_ch0]);
        else
            temp_ch = str2double(temp_ch0);
        end
        if old_ch == temp_ch, new_ch = 0; else new_ch = 1; end
        old_ch = temp_ch;
        
        plxUnits{ii,1} = dataFile{ii};
        plxUnits{ii,2}{jj,1} = temp_ch;
        plxUnits{ii,2}{jj,2} = fn{jj}(end);
        
        channel = [channel;plxUnits{ii,2}(jj,1)];
        unit = [unit;{[num2str(plxUnits{ii,2}{jj,1}),plxUnits{ii,2}{jj,2}]}];
        
        if useLogFiles
            for cc = 1:tmp_nChamber
                if any(ismember(subHeader_var.channels{cc},temp_ch))
                    tmp_area = [tmp_area;subHeader_var.Properties.RowNames(cc)];
                    tmp_nUnit(cc) = tmp_nUnit(cc)+1;
                    if new_ch, tmp_nChannel(cc) = tmp_nChannel(cc) + 1; end
                end
            end
        else
            switch tmp_nChamber
                case 1
                    for cc = 1:2
                        if any(ismember(df_electrodeInChambers(cc,:),temp_ch))
                            tmp_area = [tmp_area;df_area(cc)];
                            tmp_nUnit = tmp_nUnit+1;
                            if new_ch, tmp_nChannel = tmp_nChannel + 1; end
                        end
                    end
                case 2
                    for cc = 1:tmp_nChamber
                        if any(ismember(df_electrodeInChambers(cc,:),temp_ch))
                            tmp_area = [tmp_area;df_area(cc)];
                            tmp_nUnit(cc) = tmp_nUnit(cc)+1;
                            if new_ch, tmp_nChannel(cc) = tmp_nChannel(cc) + 1; end
                        end
                    end
            end
        end
        tmp_lfp = {[]};
        LFPch = [LFPch;tmp_lfp];
    end
    area_3 = [area_3;tmp_area];
    % 2nd level data from single unit data
    nChannel = [nChannel;num2cell(tmp_nChannel)];
    nUnit = [nUnit;num2cell(tmp_nUnit)];
    if ~useLogFiles % check again
        switch tmp_nChamber
            case 1, area = [area;unique(tmp_area)];
            case 2, area = [area;df_area];
        end
    end
    
    plxUnitsSum{ii}(:,1) = unique(cell2mat(plxUnits{ii,2}(:,1)));
    for kk = 1:length(plxUnitsSum{ii})
        plxUnitsSum{ii}(kk,2) = length(find(plxUnitsSum{ii}(kk,1) == cell2mat(plxUnits{ii,2}(:,1))));
    end
end

% save plxUnits
save([path 'plxUnits'],'plxUnits')

empty_field = LFPch;
    
% write xls file (sheet1: sessions, sheet2: chambers, sheet3: units)
new_table = cell(1,nSheet);
fin_table = cell(1,nSheet);
new_table{1} = cell2table([date,rec,dataFile,typeOfFile,monkey,task,task_variation,general_comment],'VariableNames',df_header_dataFiles);
% table2 = cell2table([date_2,dataFile_2,area,drive,head,Xpos,Ypos,ZDepth,guideTube,nElectrode,nChannel,nUnit],'VariableNames',df_header_chambers);
% table3 = cell2table([date_3,rec_3,dataFile_3,monkey_3,area_3,task_3,task_variation_3,channel,unit,LFPch,empty_field,empty_field,empty_field,empty_field,empty_field,empty_field,empty_field,empty_field],'VariableNames',df_header_units);

new_table{2} = cell2table([date_2,dataFile_2,area,drive,head,matrixCoor,ZDepth,guideTube,nElectrode,nChannel,nUnit,comment],'VariableNames',df_header_chambers);
new_table{3} = cell2table([date_3,rec_3,dataFile_3,monkey_3,area_3,task_3,task_variation_3,channel,unit,LFPch,empty_field,empty_field,empty_field,empty_field,empty_field,empty_field,empty_field,empty_field,empty_field],'VariableNames',df_header_units);

if ~options.overwrite
    % append new data files to existing table
    new_datafile = setdiff(new_table{1}.DataFile,old_table{1}.DataFile);
    for sh = 1:nSheet
        new_table{sh}.DataFile = categorical(new_table{sh}.DataFile);
        fin_table{sh} = old_table{sh};
        for nn = 1:length(new_datafile)
            fin_table{sh} = [fin_table{sh};new_table{sh}(new_table{sh}.DataFile == new_datafile(nn),:)];
        end
    end
    fin_table{3} = sortrows(fin_table{3},{'DataFile','Unit'});
    fin_table{2} = sortrows(fin_table{2},{'DataFile','Area'});
    fin_table{1} = sortrows(fin_table{1},'DataFile');
else
    for sh = 1:nSheet
        fin_table{sh} = new_table{sh};
    end
end

writetable(fin_table{3},recordingFullFile,'Sheet','units');
writetable(fin_table{2},recordingFullFile,'Sheet','chambers');
writetable(fin_table{1},recordingFullFile,'Sheet','dataFiles');

hExcel = actxserver('Excel.Application');
hWB = hExcel.Workbooks.Open(recordingFullFile); 

% remove default sheets
sheetName = 'Sheet'; % EN: Sheet, DE: Tabelle, etc. (Lang. dependent)
try
      % Throws an error if the sheets do not exist.
      hWB.Worksheets.Item([sheetName '1']).Delete;
      hWB.Worksheets.Item([sheetName '2']).Delete;
      hWB.Worksheets.Item([sheetName '3']).Delete;
catch
      ; % Do nothing.
end

% set auto-width
hWB.Worksheets.Item('units').Cells.EntireColumn.AutoFit;
hWB.Worksheets.Item('chambers').Cells.EntireColumn.AutoFit;
hWB.Worksheets.Item('dataFiles').Cells.EntireColumn.AutoFit;

% set colors in chamber & unit sheets
alphabets = ('A':'Z').';

% header
color_sheets = {'dataFiles','4F4F2F';'chambers','8B704A';'units','14698B'}; 
nCol = [size(fin_table{1},2);size(fin_table{2},2);size(fin_table{3},2)]; 
for ss = 1:size(color_sheets,1)
    hWB.Worksheets.Item(color_sheets{ss,1}).Range(['A1:' alphabets(nCol(ss)) '1']).HorizontalAlignment = 3;
    hWB.Worksheets.Item(color_sheets{ss,1}).Range(['A1:' alphabets(nCol(ss)) '1']).Interior.Color = hex2dec(color_sheets{ss,2});
    hWB.Worksheets.Item(color_sheets{ss,1}).Range(['A1:' alphabets(nCol(ss)) '1']).font.Color = hex2dec('FFFFFF');
    hWB.Worksheets.Item(color_sheets{ss,1}).Range(['A1:' alphabets(nCol(ss)) '1']).font.Bold = 1;
end

% chamber sheet
% color_areas = {'PMd','C1FFC1';'PRR','7AA0FF'};
color_areas = {'PMd','DBDCF2';'PRR','DEF1EB'};
fin_table{2}.Area = categorical(fin_table{2}.Area);

for cc = 1:size(color_areas,2)
    table2_ind = find(fin_table{2}.Area == color_areas{cc,1})';
    table2_range = [];
    for ii = table2_ind
        table2_range = [alphabets(1) num2str(ii+1) ':' alphabets(nCol(2)) num2str(ii+1)];
        hWB.Worksheets.Item(color_sheets{2,1}).Range(table2_range).Interior.Color = hex2dec(color_areas{cc,2});
    end
end

% unit sheet
fin_table{3}.Area = categorical(fin_table{3}.Area);
for cc = 1:size(color_areas,2)
    table3_ind = find(fin_table{3}.Area == color_areas{cc,1})';
    table3_range = [];
    last_ii = 0;
    for ii = table3_ind
        if last_ii && ii == (last_ii+1)
            length_last_ii = length(num2str(ii))+1;
            table3_range = [table3_range(1:end-length_last_ii) [alphabets(nCol(3)) num2str(ii+1)]];
            if ii == table3_ind(end), hWB.Worksheets.Item(color_sheets{3,1}).Range(table3_range).Interior.Color = hex2dec(color_areas{cc,2}); end
        else
            if last_ii > 0, hWB.Worksheets.Item(color_sheets{3,1}).Range(table3_range).Interior.Color = hex2dec(color_areas{cc,2}); end
            table3_range = [alphabets(1) num2str(ii+1) ':' alphabets(nCol(3)) num2str(ii+1)];
        end
        last_ii = ii;
    end
end

hWB.Save();
hExcel.Quit();
hExcel.delete();
