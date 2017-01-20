function new_filename = getSpkTrain_hdf52mat(data_folder,mat_folder,sel_monkey,pool_monkey,sel_area,sel_taskvar,pool_taskvar,align_mode)

% get_spkTrain.m gets spk data from HDF5 files in data_folder 
% according to the selected monkeys, areas (sel_monkey/sel_area), 
% alignments (align_mode), and tuning properties (selTuning) then  
% save 'spk_table' (matlab data table format) for each monkey/area/alignment 
% separately in 'mat_folder' (generated inside data_folder)

delta_t = 0.01;

%%%%%%%%%%%%%%%%%%%%%
% alignment presets %
%%%%%%%%%%%%%%%%%%%%%
% cue space
options.lateMem_cue.name  = 'lateMem_cue';
options.lateMem_cue.alignMode = 'lateMem300';
options.lateMem_cue.dtun  = 'CueDir';           % CueDir | BiasDir | MoveDir
options.lateMem_cue.filter.Success  = 1;        % selected trials for tuning (can be any task variables)
options.lateMem_cue.sorter.BiasType = [0 1 2 3 4 5 6];  % sort trial types into different conditions (can be any task variables)

options.lateMem_cue_cf.name  = 'lateMem_cue_cf';
options.lateMem_cue_cf.alignMode = 'lateMem300';
options.lateMem_cue_cf.dtun  = 'CueDir';           % CueDir | BiasDir | MoveDir
options.lateMem_cue_cf.filter.Success  = 1;        % selected trials for tuning (can be any task variables)
options.lateMem_cue_cf.sorter.BiasType = [0 1 2 3 4 5 6];  % sort trial types into different conditions (can be any task variables)
options.lateMem_cue_cf.sorter.CuingType   = [1 2];
options.lateMem_cue_cf.sorter.Follow      = [0 1]; 

% motor goal space 
options.lateMem_bias.name  = 'lateMem_bias';
options.lateMem_bias.alignMode = 'lateMem300';
options.lateMem_bias.dtun  = 'BiasDir';           % CueDir | BiasDir | MoveDir
options.lateMem_bias.filter.Success  = 1;        % selected trials for tuning (can be any task variables)
options.lateMem_bias.sorter.BiasTypeAbs = [0 1 2 3];  % sort trial types into different conditions (can be any task variables)

options.lateMem_bias_cf.name  = 'lateMem_bias_cf';
options.lateMem_bias_cf.alignMode = 'lateMem300';
options.lateMem_bias_cf.dtun  = 'BiasDir';           % CueDir | BiasDir | MoveDir
options.lateMem_bias_cf.filter.Success  = 1;        % selected trials for tuning (can be any task variables)
options.lateMem_bias_cf.sorter.BiasTypeAbs = [0 1 2 3];  % sort trial types into different conditions (can be any task variables)
options.lateMem_bias_cf.sorter.Follow      = [0 1]; 
options.lateMem_bias_cf.sorter.CuingType   = [1 2];

options.preCue.name  = 'preCue';
options.preCue.align = 'preCue1000';
options.preCue.dtun  = 'BiasDir';           % CueDir | BiasDir | MoveDir
options.preCue.filter.Success      = 1;        % selected trials for tuning (can be any task variables)
options.preCue.sorter.BiasTypeAbs = [0 1 2 3];   % sort trial types into different conditions (can be any task variables)

options.go.name  = 'go';
options.go.align = 'go700';
options.go.dtun  = 'BiasDir';           % CueDir | BiasDir | MoveDir
options.go.filter.Success      = 1;        % selected trials for tuning (can be any task variables)
options.go.sorter.BiasTypeAbs = [0 1 2 3];
options.go.sorter.CuingType    = [1 2];
% options.go.sorter.Follow  = [0 1];   % sort trial types into different conditions (can be any task variables)

options.move.name  = 'move';
options.move.align = 'move400';
options.move.dtun  = 'BiasDir';           % CueDir | BiasDir | MoveDir
options.move.filter.Success      = 1;        % selected trials for tuning (can be any task variables)
options.move.sorter.BiasTypeAbs = [0 1 2 3];
options.move.sorter.CuingType    = [1 2];
% options.move.sorter.Follow  = [0 1];   % sort trial types into different conditions (can be any task variables)


options.move_bias_IF.name  = 'move_bias_IF';
options.move_bias_IF.alignMode = 'move100';
options.move_bias_IF.dtun  = 'BiasDir';           % CueDir | BiasDir | MoveDir
options.move_bias_IF.filter.Success  = 1;        % selected trials for tuning (can be any task variables)
options.move_bias_IF.sorter.BiasTypeAbs = [0 1 2 3];  % sort trial types into different conditions (can be any task variables)
options.move_bias_IF.filter.Follow  = 1;
options.move_bias_IF.filter.CuingType  = 1;

options.move_bias_FF.name  = 'move_bias_FF';
options.move_bias_FF.alignMode = 'move100';
options.move_bias_FF.dtun  = 'BiasDir';           % CueDir | BiasDir | MoveDir
options.move_bias_FF.filter.Success  = 1;        % selected trials for tuning (can be any task variables)
options.move_bias_FF.sorter.BiasTypeAbs = [0 1 2 3];  % sort trial types into different conditions (can be any task variables)
options.move_bias_FF.filter.Follow  = 1;
options.move_bias_FF.filter.CuingType  = 2;

options.move_bias_IA.name  = 'move_bias_IA';
options.move_bias_IA.alignMode = 'move100';
options.move_bias_IA.dtun  = 'BiasDir';           % CueDir | BiasDir | MoveDir
options.move_bias_IA.filter.Success  = 1;        % selected trials for tuning (can be any task variables)
options.move_bias_IA.sorter.BiasTypeAbs = [0 1 2 3];  % sort trial types into different conditions (can be any task variables)
options.move_bias_IA.filter.Follow  = 0;
options.move_bias_IA.filter.CuingType  = 1;

options.move_bias_FA.name  = 'move_bias_FA';
options.move_bias_FA.alignMode = 'move100';
options.move_bias_FA.dtun  = 'BiasDir';           % CueDir | BiasDir | MoveDir
options.move_bias_FA.filter.Success  = 1;        % selected trials for tuning (can be any task variables)
options.move_bias_FA.sorter.BiasTypeAbs = [0 1 2 3];  % sort trial types into different conditions (can be any task variables)
options.move_bias_FA.filter.Follow  = 0;
options.move_bias_FA.filter.CuingType  = 2;

all_options = fieldnames(options);

% check whether the alignment modes are all valid & generate setttings
settings = cell(numel(align_mode),1);
for am = 1:numel(align_mode)
    if any(ismember(all_options,align_mode{am}))
        settings{am} = options.(align_mode{am});
    else
        disp([align_mode{am} 'is not a valid alignment mode'])
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load hdf5 - create spk_table - save %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for mm = 1:numel(sel_monkey)
    first_m = 1;
    for am = 1:numel(align_mode)
        first_am = 1;
        for ar = 1:numel(sel_area)
            save_filename = ['spkTrain_' sel_monkey{mm} '_BiasProb_' sel_area{ar} '_' settings{am}.name];
            if ~exist([mat_folder save_filename '.mat'],'file')
                if first_m
                    hdf_folder = [data_folder sel_monkey{mm} filesep 'hdf5' filesep];
                    [tc,spk] = load_folder(hdf_folder,'trial','.tc.hdf5','spike','.spk.hdf5');
                    td = tc.get_trial_data();
                    unit_info = spk.get_unit_info();
                    unit_table = struct2table(unit_info);
                    nUnit = size(unit_info.dataset_name,1);
                    first_m = 0;
                end
                
                if first_am
                    if ~any(ismember(fieldnames(td.variables),settings{am}.dtun)), error([settings{am}.dtun ' is not a valid tuning parameter. No matched trial variables found']); end
                    % define wdt (time window for directional tuning) and align (alignment of neural data)
                    switch settings{am}.alignMode
                        case 'preCue1000'                % 1000 ms after pre-cue onset
                            wdt    = [0, 1000];
                            align  = 'CSO';
                        case 'go700'                       % 700 ms after go cue
                            wdt    = [0, 700];
                            align  = 'C2T';
                        case 'move400'                     % 400ms before + after movement onset
                            wdt    = [-400, 400];
                            align  = 'MOV';
                        case 'lateMem300'               % 300 ms before go cue
                            wdt    = [-300, 0];
                            align  = 'C2T';
                        case 'lateMem500'               % 500 ms before go cue
                            wdt    = [-500, 0];
                            align  = 'C2T';
                        otherwise, error('No valid mode specified!')
                    end
                    
                    if strcmp(align,'MOV')
                        td.times.MOV = td.times.C2T + td.variables.ReaTime;
                        td.dates.MOV = td.dates.C2T + ms2d(td.variables.ReaTime);
                    end
                    
                    ref = td.dates.(align);
                    w_start = ref + ms2d(wdt(1));
                    w_end   = ref + ms2d(wdt(2));
                    td.variables.win_sec = round(d2ms(w_end - w_start))/1000; % time window length in second
                    % convert to data tables
                    trial_table = struct2table(td.variables);
                    first_am = 0;
                    time_bin = wdt(1):wdt(2);
                end
                %%%%%%%%%%%%%%%%%%
                % get spike data %
                %%%%%%%%%%%%%%%%%%
                selected_units = strcmp(unit_info.Area,sel_area{ar}) & strcmp(unit_info.TaskVariation,sel_taskvar{tt}) & unit_info.tuning_lateMem2 & unit_info.valide;
                dataset_names = unit_info.dataset_name(selected_units);
                [spikes,unit_ID,trial_ID] = spk.get_all_spikes(selected_units,w_start,w_end);
                
                trial_table_tmp = trial_table(trial_ID,:);
                unit_table_tmp = unit_table(unit_ID,:);
                
                % make spike bins
                int_time_bin = 1:(length(time_bin)-1) * (0.001/delta_t);
                spk_train = cellfun(@(s)histc(s,int_time_bin),spikes,'UniformOutput',0);
                
                %%%%%%%%%%%%%%%%%%%%%%%
                % filtering & sorting %
                %%%%%%%%%%%%%%%%%%%%%%%
                
                trialVar_names = trial_table_tmp.Properties.VariableNames;
                ntrial = height(trial_table_tmp);
                
                [filter,filter_names,nfilter] = filter_trials(trial_table_tmp,settings{am}.filter);
                
                if ~isempty(settings{am}.sorter)
                    [sorter_table,sorter_names,nsorter,sorter_Varname,total_cond] = sorter_trials(trial_table_tmp,settings{am}.sorter);
                else
                    sorter_names = {'all'};
                    total_cond = 1;
                    sorter_table = ones(height(trial_table_tmp),total_cond);
                    sorter_Varname = {'all'};
                    nsorter = 0;
                end
                sorter_table = array2table(sorter_table,'VariableNames',sorter_Varname);
                
                header = [{'unitName','prefDir','maxAct','maxDir','spikeDensity','cueDir','biasDir','moveDir','fileIndex','trialIndex','condName'}, sorter_names];
                spkTrain_table = cell2table(cell(sum(filter),length(header)),'VariableNames',header);
                spkTrain_table.unitName = unit_table_tmp.dataset_name(filter);
                
                spkTrain_table.prefDir = unit_table_tmp.tuning_lateMem2_angle(filter);
                spkTrain_table.maxAct = unit_table_tmp.tuning_lateMem2_maxAct(filter);
                spkTrain_table.maxDir = unit_table_tmp.tuning_lateMem2_maxDir(filter);
                
                spkTrain_table.spikeTrain = spk_train(filter);
                spkTrain_table.spikeTrain = cell2mat(spkTrain_table.spikeTrain);
                
                spkTrain_table.cueDir = trial_table_tmp.CueDir(filter);
                spkTrain_table.biasDir = trial_table_tmp.BiasDir(filter);
                spkTrain_table.moveDir = trial_table_tmp.HDirMv(filter);
                spkTrain_table.trialIndex = trial_table_tmp.trial_index(filter);
                spkTrain_table.fileIndex = trial_table_tmp.file_index(filter);
                
                for ss = 1:length(sorter_names)
                    spkTrain_table.(sorter_names{ss}) = trial_table_tmp.(sorter_names{ss})(filter);
                end
                
                spkTrain_table.condName = spkTrain_table.unitName;
                for cc = 1:total_cond
                    spkTrain_table.condName(logical(sorter_table{:,cc}(filter))) = repmat(sorter_Varname(cc),sum(sorter_table{:,cc}(filter)),1);
                end
                
                setting = settings{am};
                
                disp(['Saving spk train data to: ',mat_folder save_filename '.mat']);
                save([mat_folder save_filename  '.mat'], 'spkTrain_table' ,'sel_unit_table','setting');
                clear spkTrain_table
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%
% pooling data & save %
%%%%%%%%%%%%%%%%%%%%%%%
new_filename = [];
for am = 1:numel(align_mode)
    for ar = 1:numel(sel_area)
        if pool_monkey
            tmp_spkTrain_table = [];
            tmp_filename = ['spkTrain_allMonkeys_BiasProb_' sel_area{ar} '_' settings{am}.name];
            if ~exist([mat_folder tmp_filename '.mat'],'file')
                last_f = 0; last_t = 0;
                for mm = 1:numel(sel_monkey)
                    save_filename = ['spkTrain_' sel_monkey{mm} '_BiasProb_' sel_area{ar} '_' settings{am}.name];
                    disp(['Loading spk train data from: ',mat_folder save_filename '.mat']);
                    load([mat_folder save_filename '.mat']);
                    spkTrain_table.fileIndex = spkTrain_table.fileIndex + last_f;
                    spkTrain_table.trialIndex = spkTrain_table.trialIndex + last_t;
                    tmp_spkTrain_table = [tmp_spkTrain_table;spkTrain_table];
                    tmp_spkTrain_table = sortrows(tmp_spkTrain_table,{'fileIndex','trialIndex'});
                    last_f = tmp_spkTrain_table.fileIndex(end);
                    last_t = tmp_spkTrain_table.trialIndex(end);
                    clear spkTrain_table
                end
                spkTrain_table = tmp_spkTrain_table;
                disp(['Saving pooled spk train data to: ',mat_folder tmp_filename '.mat']);
                save([mat_folder tmp_filename '.mat'], 'spkTrain_table',...
                    'time_bin','align','sorter_Varname','setting','-v7.3');
            end
            new_filename = [new_filename; cellstr(tmp_filename)];
        else
            tmp_filename = cell(numel(sel_monkey),1);
            for mm = 1:numel(sel_monkey)
                save_filename = ['spkTrain_' sel_monkey{mm} '_BiasProb_' sel_area{ar} '_' settings{am}.name];
                tmp_filename{mm} = save_filename;
            end
            new_filename = [new_filename; tmp_filename];
        end
    end
end