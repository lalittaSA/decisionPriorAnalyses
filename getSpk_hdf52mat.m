function new_filename = getSpk_hdf52mat(data_folder,mat_folder,sel_monkey,pool_monkey,sel_area,sel_taskVar,pool_taskVar,align_mode)

% get_spk.m gets spk data from HDF5 files in data_folder 
% according to the selected monkeys, areas (sel_monkey/sel_area), 
% alignments (align_mode), and tuning properties (selTuning) then  
% save 'spk_table' (matlab data table format) for each monkey/area/alignment 
% separately in 'mat_folder' (generated inside data_folder)

%%%%%%%%%%%%%%%%%%%%%
% alignment presets %
%%%%%%%%%%%%%%%%%%%%%
% cue space
options.lateMem_cue.name  = 'lateMem_cue';
options.lateMem_cue.alignMode = 'lateMem300';
options.lateMem_cue.dtun  = 'cueDir';           % CueDir | BiasDir | MoveDir
options.lateMem_cue.filter.Success  = 1;        % selected trials for tuning (can be any task variables)
options.lateMem_cue.sorter.BiasType = [0 1 2 3 4 5 6];  % sort trial types into different conditions (can be any task variables)

options.lateMem_cue_cf.name  = 'lateMem_cue_cf';
options.lateMem_cue_cf.alignMode = 'lateMem300';
options.lateMem_cue_cf.dtun  = 'cueDir';           % CueDir | BiasDir | MoveDir
options.lateMem_cue_cf.filter.Success  = 1;        % selected trials for tuning (can be any task variables)
options.lateMem_cue_cf.sorter.BiasType = [0 1 2 3 4 5 6];  % sort trial types into different conditions (can be any task variables)
options.lateMem_cue_cf.sorter.CuingType   = [1 2];
options.lateMem_cue_cf.sorter.Follow      = [0 1]; 

% motor goal space 
options.lateMem_bias.name  = 'lateMem_bias';
options.lateMem_bias.alignMode = 'lateMem300';
options.lateMem_bias.dtun  = 'biasDir';           % CueDir | BiasDir | MoveDir
options.lateMem_bias.filter.Success  = 1;        % selected trials for tuning (can be any task variables)
options.lateMem_bias.sorter.BiasTypeAbs = [0 1 2 3];  % sort trial types into different conditions (can be any task variables)

options.lateMem_bias_cf.name  = 'lateMem_bias_cf';
options.lateMem_bias_cf.alignMode = 'lateMem300';
options.lateMem_bias_cf.dtun  = 'biasDir';           % CueDir | BiasDir | MoveDir
options.lateMem_bias_cf.filter.Success  = 1;        % selected trials for tuning (can be any task variables)
options.lateMem_bias_cf.sorter.BiasTypeAbs = [0 1 2 3];  % sort trial types into different conditions (can be any task variables)
options.lateMem_bias_cf.sorter.Follow      = [0 1]; 
options.lateMem_bias_cf.sorter.CuingType   = [1 2];

options.move_bias_IF.name  = 'move_bias_IF';
options.move_bias_IF.alignMode = 'move100';
options.move_bias_IF.dtun  = 'biasDir';           % CueDir | BiasDir | MoveDir
options.move_bias_IF.filter.Success  = 1;        % selected trials for tuning (can be any task variables)
options.move_bias_IF.sorter.BiasTypeAbs = [0 1 2 3];  % sort trial types into different conditions (can be any task variables)
options.move_bias_IF.filter.Follow  = 1;
options.move_bias_IF.filter.CuingType  = 1;

options.move_bias_FF.name  = 'move_bias_FF';
options.move_bias_FF.alignMode = 'move100';
options.move_bias_FF.dtun  = 'biasDir';           % CueDir | BiasDir | MoveDir
options.move_bias_FF.filter.Success  = 1;        % selected trials for tuning (can be any task variables)
options.move_bias_FF.sorter.BiasTypeAbs = [0 1 2 3];  % sort trial types into different conditions (can be any task variables)
options.move_bias_FF.filter.Follow  = 1;
options.move_bias_FF.filter.CuingType  = 2;

options.move_bias_IA.name  = 'move_bias_IA';
options.move_bias_IA.alignMode = 'move100';
options.move_bias_IA.dtun  = 'biasDir';           % CueDir | BiasDir | MoveDir
options.move_bias_IA.filter.Success  = 1;        % selected trials for tuning (can be any task variables)
options.move_bias_IA.sorter.BiasTypeAbs = [0 1 2 3];  % sort trial types into different conditions (can be any task variables)
options.move_bias_IA.filter.Follow  = 0;
options.move_bias_IA.filter.CuingType  = 1;

options.move_bias_FA.name  = 'move_bias_FA';
options.move_bias_FA.alignMode = 'move100';
options.move_bias_FA.dtun  = 'biasDir';           % CueDir | BiasDir | MoveDir
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
            for tt = 1:numel(sel_taskVar)
                save_filename = ['spk_' sel_monkey{mm} '_' sel_taskVar{tt} '_' sel_area{ar} '_' settings{am}.name '_tuned'];
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
                        % define wdt (time window for directional tuning) and align (alignment of neural data)
                        switch settings{am}.alignMode
                            case 'move100'               % 100 ms around movement initiation
                                wdt    = [-50, 50];
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
                    end
                    %%%%%%%%%%%%%%%%%%
                    % get spike data %
                    %%%%%%%%%%%%%%%%%%
                    selected_units = strcmp(unit_info.Area,sel_area{ar}) & strcmp(unit_info.TaskVariation,sel_taskVar{tt}) & unit_info.tuning_lateMem2 & unit_info.valide;
                    dataset_names = unit_info.dataset_name(selected_units);
                    [spikes,unit_ID,trial_ID] = spk.get_all_spikes(selected_units,w_start,w_end);
                    
                    trial_table_tmp = trial_table(trial_ID,:);
                    unit_table_tmp = unit_table(unit_ID,:);
                    
                    %Compute n_spikes & firing rates for tuning stat
                    nspikes = cellfun(@(s)numel(s),spikes);
                    firing = nspikes./trial_table_tmp.win_sec;
                    
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
                    
                    spk_table = [];
                    used_unit = zeros(total_cond,1);
                    used_ind = false(nUnit,1);
                    for uu = 1:nUnit
                        unit_filt = filter & ismember(unit_table_tmp.dataset_name,unit_info.dataset_name(uu));
                        for cc = 1:total_cond
                            unit_sel = unit_filt & sorter_table{:,cc};
                            if sum(unit_sel)
                                used_ind(uu) = 1;
                                used_unit(cc) = used_unit(cc) + 1;
                                tmp_unitName = unit_table_tmp.dataset_name(unit_sel);
                                tmp_prefDir = num2cell(unit_table_tmp.tuning_lateMem2_angle(unit_sel));
                                tmp_maxAct = num2cell(unit_table_tmp.tuning_lateMem2_maxAct(unit_sel));
                                tmp_maxDir = num2cell(unit_table_tmp.tuning_lateMem2_maxDir(unit_sel));
                                tmp_firing = num2cell(firing(unit_sel));
                                tmp_cueDir = num2cell(trial_table_tmp.CueDir(unit_sel));
                                tmp_biasDir = num2cell(trial_table_tmp.BiasDir(unit_sel));
                                tmp_moveDir = num2cell(trial_table_tmp.HDirMv(unit_sel));
                                tmp_trialID = num2cell(trial_table_tmp.trial_index(unit_sel));
                                tmp_fileID = num2cell(trial_table_tmp.file_index(unit_sel));
                                tmp_VarName = repmat(sorter_Varname(cc),sum(unit_sel),1);
                                tmp_cond = [];
                                for st = 1:nsorter
                                    tmp_cond = [tmp_cond,num2cell(trial_table_tmp.(sorter_names{st})(unit_sel))];
                                end
                                cond_table = cell2table(tmp_cond,'VariableNames',sorter_names);
                                tmp_spk_table = [cell2table([tmp_unitName,tmp_prefDir,tmp_maxAct,tmp_maxDir,tmp_firing,tmp_cueDir,tmp_biasDir,tmp_moveDir,tmp_fileID,tmp_trialID,tmp_VarName],...
                                    'VariableNames',{'unitName','prefDir','maxAct','maxDir','firingRate','cueDir','biasDir','moveDir','fileIndex','trialIndex','condName'}),...
                                    cond_table];
                                spk_table = [spk_table; tmp_spk_table];
                            end
                        end
                    end
                    sel_unit_table = unit_table(used_ind,:);

                    setting = settings{am};
                    disp(['Saving spk data to: ',mat_folder save_filename '.mat']);
                    save([mat_folder save_filename  '.mat'], 'spk_table' ,'sel_unit_table','setting');
                    clear spk_table
                end
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
        if pool_monkey && pool_taskVar
            tmp_spk_table = [];
            tmp_sel_unit_table = [];
            tmp_filename = ['spk_allMonkeys_allTaskVar_' sel_area{ar} '_' settings{am}.name '_tuned'];
            if ~exist([mat_folder tmp_filename '.mat'],'file')
                last_f = 0; last_t = 0;
                for mm = 1:numel(sel_monkey)
                    for tt = 1:numel(sel_taskVar)
                        save_filename = ['spk_' sel_monkey{mm} '_' sel_taskVar{tt} '_' sel_area{ar} '_' settings{am}.name '_tuned'];
                        disp(['Loading spk data from: ',mat_folder save_filename '.mat']);
                        load([mat_folder save_filename '.mat']);
                        spk_table.fileIndex = spk_table.fileIndex + last_f;
                        spk_table.trialIndex = spk_table.trialIndex + last_t;
                        tmp_spk_table = [tmp_spk_table; spk_table];
                        tmp_sel_unit_table = [tmp_sel_unit_table;sel_unit_table];
                    end
                    tmp_spk_table = sortrows(tmp_spk_table,{'fileIndex','trialIndex'});
                    last_f = tmp_spk_table.fileIndex(end);
                    last_t = tmp_spk_table.trialIndex(end);
                    clear spk_table sel_unit_table
                end
                spk_table = tmp_spk_table;
                sel_unit_table = tmp_sel_unit_table;
                disp(['Saving pooled spk data to: ',mat_folder tmp_filename '.mat']);
                save([mat_folder tmp_filename  '.mat'], 'spk_table', 'sel_unit_table','setting');
            end
            new_filename = [new_filename; cellstr(tmp_filename)];
        elseif ~pool_monkey && pool_taskVar
            tmp_filename = cell(numel(sel_monkey),1);
            for mm = 1:numel(sel_monkey)
                tmp_filename{mm} = ['spk_' sel_monkey{mm} '_allTaskVar_' sel_area{ar} '_' settings{am}.name '_tuned'];
                if ~exist([mat_folder tmp_filename{mm} '.mat'],'file')
                    tmp_spk_table = [];
                    tmp_sel_unit_table = [];
                    for tt = 1:numel(sel_taskVar)
                        last_f = 0; last_t = 0;
                        save_filename = ['spk_' sel_monkey{mm} '_' sel_taskVar{tt} '_' sel_area{ar} '_' settings{am}.name '_tuned'];
                        disp(['Loading spk data from: ',mat_folder save_filename '.mat']);
                        load([mat_folder save_filename '.mat']);
                        spk_table.fileIndex = spk_table.fileIndex + last_f;
                        spk_table.trialIndex = spk_table.trialIndex + last_t;
                        tmp_spk_table = [tmp_spk_table;spk_table];
                        tmp_sel_unit_table = [tmp_sel_unit_table;sel_unit_table];
                        tmp_spk_table = sortrows(tmp_spk_table,{'fileIndex','trialIndex'});
                        last_f = tmp_spk_table.fileIndex(end);
                        last_t = tmp_spk_table.trialIndex(end);
                        clear spk_table sel_unit_table
                    end
                    spk_table = tmp_spk_table;
                    sel_unit_table = tmp_sel_unit_table;
                    disp(['Saving pooled spk data to: ',mat_folder tmp_filename{mm} '.mat']);
                    save([mat_folder tmp_filename{mm} '.mat'], 'spk_table', 'sel_unit_table','setting');
                end
            end
            new_filename = [new_filename; tmp_filename];
        elseif pool_monkey && ~pool_taskVar
            tmp_filename = cell(numel(sel_taskVar),1);
            for tt = 1:numel(sel_taskVar)
                tmp_filename{tt} = ['spk_allMonkeys_' sel_taskVar{tt} '_' sel_area{ar} '_' settings{am}.name '_tuned'];
                if ~exist([mat_folder tmp_filename{tt} '.mat'],'file')
                    last_f = 0; last_t = 0;
                    tmp_spk_table = [];
                    tmp_sel_unit_table = [];
                    for mm = 1:numel(sel_monkey)
                        save_filename = ['spk_' sel_monkey{mm} '_' sel_taskVar{tt} '_' sel_area{ar} '_' settings{am}.name '_tuned'];
                        disp(['Loading spk data from: ',mat_folder save_filename '.mat']);
                        load([mat_folder save_filename '.mat']);
                        spk_table.fileIndex = spk_table.fileIndex + last_f;
                        spk_table.trialIndex = spk_table.trialIndex + last_t;
                        tmp_spk_table = [tmp_spk_table;spk_table];
                        tmp_sel_unit_table = [tmp_sel_unit_table;sel_unit_table];
                        tmp_spk_table = sortrows(tmp_spk_table,{'fileIndex','trialIndex'});
                        last_f = tmp_spk_table.fileIndex(end);
                        last_t = tmp_spk_table.trialIndex(end);
                        clear spk_table sel_unit_table
                    end
                    spk_table = tmp_spk_table;
                    sel_unit_table = tmp_sel_unit_table;
                    disp(['Saving pooled spk data to: ',mat_folder tmp_filename{tt} '.mat']);
                    save([mat_folder tmp_filename{tt} '.mat'], 'spk_table', 'sel_unit_table','setting');
                end
            end
            new_filename = [new_filename; tmp_filename];
        else
            tmp_filename = cell(numel(sel_monkey)*numel(sel_taskVar),1);
            mt = 1;
            for mm = 1:numel(sel_monkey)
                for tt = 1:numel(sel_taskVar)
                     save_filename = ['spk_' sel_monkey{mm} '_' sel_taskVar{tt} '_' sel_area{ar} '_' settings{am}.name '_tuned'];
                     tmp_filename{mt} = save_filename;
                     mt = mt+1;
                end
            end
            new_filename = [new_filename; tmp_filename];
        end
    end
end