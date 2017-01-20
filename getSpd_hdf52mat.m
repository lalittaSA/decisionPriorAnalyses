function new_filename = getSpd_hdf52mat(data_folder,mat_folder,pool_monkey,sel_monkey,sel_area,align_mode,spd_type,sel_tuning)

% get_spd.m gets spd data (spd_type) from HDF5 files in data_folder) 
% according to the selected monkeys, areas (selMonkey/selArea), 
% alignments (align_mode), and tuning properties (selTuning) then  
% save 'spd_table' (matlab data table format) for each monkey/area/alignment 
% separately in 'mat_folder' (generated inside data_folder)

%%%%%%%%%%%%%%%%%%%%%
% alignment presets %
%%%%%%%%%%%%%%%%%%%%%
options.preCue.name  = 'preCue';
options.preCue.align = 'preCue1000';
options.preCue.dtun  = 'biasDir';           % CueDir | BiasDir | MoveDir
options.preCue.filter.Success      = 1;        % selected trials for tuning (can be any task variables)
options.preCue.sorter.BiasTypeAbs = [0 1 2 3];   % sort trial types into different conditions (can be any task variables)

% test Hump
options.preCueSacc.name  = 'preCueSacc';
options.preCueSacc.align = 'sacc2F1000';
options.preCueSacc.dtun  = 'biasDir';           % CueDir | BiasDir | MoveDir
options.preCueSacc.filter.Success      = 1;        % selected trials for tuning (can be any task variables)
options.preCueSacc.sorter.BiasTypeAbs = [0 1 2 3];   % sort trial types into different conditions (can be any task variables)

options.preCue_cue.name  = 'preCue_cue';
options.preCue_cue.align = 'preCue1000';
options.preCue_cue.dtun  = 'cueDir';           % CueDir | BiasDir | MoveDir
options.preCue_cue.filter.Success      = 1;        % selected trials for tuning (can be any task variables)
options.preCue_cue.sorter.BiasType = [0 1 2 3 4 5 6];   % sort trial types into different conditions (can be any task variables)

options.preCue_cue_cf.name  = 'preCue_cue_cf';
options.preCue_cue_cf.align = 'preCue1000';
options.preCue_cue_cf.dtun  = 'cueDir';           % CueDir | BiasDir | MoveDir
options.preCue_cue_cf.filter.Success      = 1;        % selected trials for tuning (can be any task variables)
options.preCue_cue_cf.sorter.BiasType = [0 1 2 3 4 5 6];   % sort trial types into different conditions (can be any task variables)
options.preCue_cue_cf.sorter.CuingType    = [1 2];
options.preCue_cue_cf.sorter.Follow  = [0 1];   % sort trial types into different conditions (can be any task variables)

options.preCue_cf.name  = 'preCue_cf';
options.preCue_cf.align = 'preCue1000';
options.preCue_cf.dtun  = 'biasDir';           % CueDir | BiasDir | MoveDir
options.preCue_cf.filter.Success      = 1;        % selected trials for tuning (can be any task variables)
options.preCue_cf.sorter.BiasTypeAbs = [0 1 2 3];   % sort trial types into different conditions (can be any task variables)
options.preCue_cf.sorter.CuingType    = [1 2];
options.preCue_cf.sorter.Follow  = [0 1];   % sort trial types into different conditions (can be any task variables)

options.lateMem.name  = 'lateMem';
options.lateMem.align = 'lateMem500';
options.lateMem.dtun  = 'biasDir';           % CueDir | BiasDir | MoveDir
options.lateMem.filter.Success      = 1;        % selected trials for tuning (can be any task variables)
options.lateMem.sorter.BiasTypeAbs = [0 1 2 3];
% options.lateMem.sorter.MappingType  = [1 2];   % sort trial types into different conditions (can be any task variables)

options.lateMem_cf.name  = 'lateMem_cf';
options.lateMem_cf.align = 'lateMem500';
options.lateMem_cf.dtun  = 'biasDir';           % CueDir | BiasDir | MoveDir
options.lateMem_cf.filter.Success      = 1;        % selected trials for tuning (can be any task variables)
options.lateMem_cf.sorter.BiasTypeAbs = [0 1 2 3];
options.lateMem_cf.sorter.CuingType    = [1 2];
options.lateMem_cf.sorter.Follow  = [0 1]; 

options.go.name  = 'go';
options.go.align = 'go700';
options.go.dtun  = 'biasDir';           % CueDir | BiasDir | MoveDir
options.go.filter.Success      = 1;        % selected trials for tuning (can be any task variables)
options.go.sorter.BiasTypeAbs = [0 1 2 3];
options.go.sorter.CuingType    = [1 2];
options.go.sorter.Follow  = [0 1];   % sort trial types into different conditions (can be any task variables)

options.go_cue.name  = 'go_cue';
options.go_cue.align = 'go700';
options.go_cue.dtun  = 'cueDir';           % CueDir | BiasDir | MoveDir
options.go_cue.filter.Success      = 1;        % selected trials for tuning (can be any task variables)
options.go_cue.sorter.BiasType = [0 1 2 3 4 5 6];
options.go_cue.sorter.CuingType    = [1 2];
options.go_cue.sorter.Follow  = [0 1];   % sort trial types into different conditions (can be any task variables)

options.move.name  = 'move';
options.move.align = 'move400';
options.move.dtun  = 'biasDir';           % CueDir | BiasDir | MoveDir
options.move.filter.Success      = 1;        % selected trials for tuning (can be any task variables)
options.move.sorter.BiasTypeAbs = [0 1 2 3];
options.move.sorter.CuingType    = [1 2];
options.move.sorter.Follow  = [0 1];   % sort trial types into different conditions (can be any task variables)

options.move_cue.name  = 'move_cue';
options.move_cue.align = 'move400';
options.move_cue.dtun  = 'cueDir';           % CueDir | BiasDir | MoveDir
options.move_cue.filter.Success      = 1;        % selected trials for tuning (can be any task variables)
options.move_cue.sorter.BiasType = [0 1 2 3 4 5 6];
options.move_cue.sorter.CuingType    = [1 2];
options.move_cue.sorter.Follow  = [0 1];   % sort trial types into different conditions (can be any task variables)

options.touch.name  = 'touch';
options.touch.align = 'preTouch500';
options.touch.dtun  = 'biasDir';           % CueDir | BiasDir | MoveDir
options.touch.filter.Success      = 1;        % selected trials for tuning (can be any task variables)
options.touch.sorter.BiasTypeAbs = [0 1 2 3];
options.touch.sorter.CuingType    = [1 2];
options.touch.sorter.Follow  = [0 1];   % sort trial types into different conditions (can be any task variables)

% get all available options
all_options = fieldnames(options);

% check whether the alignment modes are all valid & generate setttings
settings = cell(numel(align_mode),1);
for am = 1:numel(align_mode)
    if any(ismember(all_options,align_mode{am}))
        settings{am} = options.(align_mode{am});
    else
        disp([align_mode{am} ' is not a valid alignment mode'])
    end
end

%%
switch sel_tuning
    case 1, file_suff = '_tuned';
    case 2, file_suff = '_tuned2';
    case 3, file_suff = '_tuned3';
    otherwise, file_suff = '';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load hdf5 - create spd_table - save %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for mm = 1:length(sel_monkey)
    first_m = 1;
    for am = 1:numel(align_mode)
        first_am = 1;
        for ar = 1:length(sel_area)          
            save_filename = ['spd' spd_type '_' sel_monkey{mm} '_allTaskVar_' sel_area{ar} '_' settings{am}.name file_suff];
            load_data   = 0;
            if exist([mat_folder save_filename '.mat'],'file')
                load_data   = 1;
            end
            if ~load_data
                tic
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
                    switch settings{am}.align
                        case 'preCue1000'                % 1000 ms after pre-cue onset
                            wdt    = [-100, 1000];
                            align  = 'CSO';
                        case 'sacc1N1000'                % 1000 ms after saccade-out onset (from fixation point)
                            wdt    = [-100, 1000];
                            align  = 'S1N';
                        case 'sacc1F1000'                % 1000 ms after saccade-out offset (from fixation point)
                            wdt    = [-100, 1000];
                            align  = 'S1F';
                        case 'sacc2N1000'                % 1000 ms after saccade-in onset (back into fixation point)
                            wdt    = [-100, 1000];
                            align  = 'S2N';
                        case 'sacc2F1000'                % 1000 ms after saccade-in offset (back into fixation point)
                            wdt    = [-100, 1000];
                            align  = 'S2F';
                        case 'lateMem500'                   % 500 ms before go cue
                            wdt    = [-500, 0];
                            align  = 'C2T';
                        case 'go700'                       % 700 ms after go cue
                            wdt    = [-100, 700];
                            align  = 'C2T';
                        case 'go1000'                       % 1000 ms after go cue
                            wdt    = [-100, 1000];
                            align  = 'C2T';
                        case 'move400'                     % 400ms before & after movement onset
                            wdt    = [-400, 400];
                            align  = 'MOV';
                        case 'move1000'                     % 1000ms after movement onset
                            wdt    = [-100, 1000];
                            align  = 'MOV';
                        case 'preTouch500'                 % 500 ms before movement offset
                            wdt    = [-500, 100];
                            align  = 'TOU';
                        case 'preTouch1000'                 % 1000 ms before movement offset
                            wdt    = [-1000, 100];
                            align  = 'TOU';
                        otherwise, error('No valid mode specified!')
                    end
                    
                    if ~iscell(align), align = cellstr(align);end
                    
                    for aa = 1:numel(align)
                        if strcmp(align{aa},'MOV')
                            td.times.MOV = td.times.C2T + td.variables.ReaTime;
                            td.dates.MOV = td.dates.C2T + ms2d(td.variables.ReaTime);
                        elseif strcmp(align{aa},'TOU')
                            td.times.TOU = td.times.C2T + td.variables.ReaTime + td.variables.MovTime;
                            td.dates.TOU = td.dates.C2T + ms2d(td.variables.ReaTime) + ms2d(td.variables.MovTime);
                        elseif strcmp(align{aa},'S1N') || strcmp(align{aa},'S1F') || strcmp(align{aa},'S2N') || strcmp(align{aa},'S2F')
                            % eye data needed to extract saccade on-/offset (test: Humphrey's data)
                            eh = load_folder(hdf_folder,'continuous','.eh.hdf5');
                            ex = eh.get_subdata_date('Eye','EyeX',td.dates.CSO,td.dates.CSF); % get eye data during cue period
                            [ey,~,~,fideye] = eh.get_subdata_date('Eye','EyeY',td.dates.CSO,td.dates.CSF);
                            tmp_screen = eh.get_attribute('Eye','screen_dist');
                            % screen distance for each trial to calculate visual angle
                            screen_dist = NaN(size(ey));
                            for tt = 1:size(ey,1)
                                if ~isnan(fideye(tt)), screen_dist(tt) = tmp_screen{fideye(tt)}; end
                            end
                            % now get the saccades
                            [SN,SF] = extract_saccades(ex,ey,screen_dist);
                            td.times.S1N = td.times.CSO + SN(:,1);
                            td.times.S1F = td.times.CSO + SF(:,1);
                            td.times.S2N = td.times.CSO + SN(:,2);
                            td.times.S2F = td.times.CSO + SF(:,2);
                            
                            td.dates.S1N = td.dates.CSO + ms2d(SN(:,1));
                            td.dates.S1F = td.dates.CSO + ms2d(SF(:,1));
                            td.dates.S2N = td.dates.CSO + ms2d(SN(:,2));
                            td.dates.S2F = td.dates.CSO + ms2d(SF(:,2));
                        end
                    end
                    
                    % get time window length for calculating firing rates
                    if numel(align) > 1
                        w_start = td.dates.(align{1});
                        w_end   = td.dates.(align{2});
                    else
                        ref = td.dates.(align{1});
                        w_start = ref + ms2d(wdt(1)-0.1); % -0.1 to get rid off numeric problem of time bins
                        w_end   = ref + ms2d(wdt(2));
                    end
                    %         win_sec = round(d2ms(w_end - w_start))/1000; % time window length in second
                    time_bin = wdt(1):wdt(2);
                    % convert to data tables
                    trial_table = struct2table(td.variables);
                    first_am = 0;
                end
                %%%%%%%%%%%%%%%%%%%%%%%
                % get spike densities %
                %%%%%%%%%%%%%%%%%%%%%%%
                switch sel_tuning
                    case 1, selected_units = strcmp(unit_info.Area,sel_area{ar}) & unit_info.valide & unit_info.tuning_lateMem2;
                    case 2, selected_units = strcmp(unit_info.Area,sel_area{ar}) & unit_info.valide & unit_info.tuning_lateMem2 & unit_info.tuning_preMove2 & (abs(unit_info.tuning_lateMem2_angle - unit_info.tuning_preMove2_angle)<90 | abs(unit_info.tuning_lateMem2_angle - unit_info.tuning_preMove2_angle)>270);
                    case 3, selected_units = strcmp(unit_info.Area,sel_area{ar}) & unit_info.valide & unit_info.tuning_lateMem2 & unit_info.tuning_preMove2 & unit_info.tuning_move2 & (abs(unit_info.tuning_lateMem2_angle - unit_info.tuning_preMove2_angle)<90 | abs(unit_info.tuning_lateMem2_angle - unit_info.tuning_preMove2_angle)>270) & (abs(unit_info.tuning_lateMem2_angle - unit_info.tuning_move2_angle)<90 | abs(unit_info.tuning_lateMem2_angle - unit_info.tuning_move2_angle)>270) & (abs(unit_info.tuning_move2_angle - unit_info.tuning_preMove2_angle)<90 | abs(unit_info.tuning_move2_angle - unit_info.tuning_preMove2_angle)>270);
                    otherwise, selected_units = strcmp(unit_info.Area,sel_area{ar}) & unit_info.valide;
                end
                dataset_names = unit_info.dataset_name(selected_units);
                %Get Continuous data
                [density,unit_ID,trial_ID,time,~] = spk.get_all_continuous(['spike_densities_' spd_type],selected_units,w_start,w_end);
                
                trial_table_tmp = trial_table(trial_ID,:);
                unit_table_tmp = unit_table(unit_ID,:);
                %align time to ref event
                %time=cellfun(@(t)t-t(1)-start_to_ref,time,'UniformOutput',false);
                time = cellfun(@(t)t-t(1)-wdt(1),time,'UniformOutput',false);
                
                
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
                spd_table = cell2table(cell(sum(filter),length(header)),'VariableNames',header);
                spd_table.unitName = unit_table_tmp.dataset_name(filter);
                
                spd_table.prefDir = unit_table_tmp.tuning_lateMem2_angle(filter);
                spd_table.maxAct = unit_table_tmp.tuning_lateMem2_maxAct(filter);
                spd_table.maxDir = unit_table_tmp.tuning_lateMem2_maxDir(filter);
                
                spd_table.spikeDensity = density(filter);
                spd_table.spikeDensity = cell2mat(spd_table.spikeDensity);
                
                spd_table.cueDir = trial_table_tmp.CueDir(filter);
                spd_table.biasDir = trial_table_tmp.BiasDir(filter);
                spd_table.moveDir = trial_table_tmp.HDirMv(filter);
                spd_table.trialIndex = trial_table_tmp.trial_index(filter);
                spd_table.fileIndex = trial_table_tmp.file_index(filter);
                
                for ss = 1:length(sorter_names)
                    spd_table.(sorter_names{ss}) = trial_table_tmp.(sorter_names{ss})(filter);
                end
                
                spd_table.condName = spd_table.unitName;
                for cc = 1:total_cond
                    spd_table.condName(logical(sorter_table{:,cc}(filter))) = repmat(sorter_Varname(cc),sum(sorter_table{:,cc}(filter)),1);
                end
                
                setting = settings{am};
                
                % clean-up: some error reaches towards the cue location in
                % choice trials survived as correct trials
                spd_table(spd_table.moveDir == spd_table.cueDir,:) = [];
                toc
                disp(['Saving spd data to: ',mat_folder save_filename '.mat']);
                save([mat_folder save_filename '.mat'], 'spd_table',...
                    'time_bin','align','sorter_Varname','setting','-v7.3');
                clear spd_table
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
            tmp_spd_table = [];
            tmp_filename = ['spd' spd_type '_allMonkeys_allTaskVar_' sel_area{ar} '_' settings{am}.name file_suff];
            if ~exist([mat_folder tmp_filename '.mat'],'file')
                last_f = 0; last_t = 0;
                for mm = 1:numel(sel_monkey)
                    save_filename = ['spd' spd_type '_' sel_monkey{mm} '_allTaskVar_' sel_area{ar} '_' settings{am}.name file_suff];
                    disp(['Loading spd data from: ',mat_folder save_filename '.mat']);
                    load([mat_folder save_filename '.mat']);
                    spd_table.fileIndex = spd_table.fileIndex + last_f;
                    spd_table.trialIndex = spd_table.trialIndex + last_t;
                    tmp_spd_table = [tmp_spd_table;spd_table];
                    tmp_spd_table = sortrows(tmp_spd_table,{'fileIndex','trialIndex'});
                    last_f = tmp_spd_table.fileIndex(end);
                    last_t = tmp_spd_table.trialIndex(end);
                    clear spd_table
                end
                spd_table = tmp_spd_table;
                disp(['Saving pooled spd data to: ',mat_folder tmp_filename '.mat']);
                save([mat_folder tmp_filename '.mat'], 'spd_table',...
                    'time_bin','align','sorter_Varname','setting','-v7.3');
            end
            new_filename = [new_filename; cellstr(tmp_filename)];
        else
            tmp_filename = cell(numel(sel_monkey),1);
            for mm = 1:numel(sel_monkey)
                save_filename = ['spd' spd_type '_' sel_monkey{mm} '_allTaskVar_' sel_area{ar} '_' settings{am}.name file_suff];
                tmp_filename{mm} = save_filename;
            end
            new_filename = [new_filename; tmp_filename];
        end
    end
end
