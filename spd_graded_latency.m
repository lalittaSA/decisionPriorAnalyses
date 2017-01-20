% spd_graded_latency.m gets spd data depending on selected time_bins / alignments then
% visualize them

%%%%%%%%%%%%%
% set paths %
%%%%%%%%%%%%%
if ismac
    data_folder = '/Volumes/Passport/dataAnalysis/data/';
    result_path = '/dataAnalysis/results/spd/';
else % in case used on Windows
    result_path = 'D:\results\hdf5\spd\';
    data_folder = 'D:\data\Decision_bias\';
end
% mat_folder = '/dataAnalysis/mat/spd/';
mat_folder = [data_folder 'mat' filesep 'spd' filesep];
figure_path = [result_path 'figures' filesep];

if ~exist(figure_path,'dir'), mkdir(figure_path); end
if ~exist(mat_folder,'dir'), mkdir(mat_folder); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select monkeys, areas, tasks %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sel_monkey = {'Hum','Kas'};
pool_monkey = 1 & numel(sel_monkey)>1;  
% 0: analyse data of all monkeys separately
% 1: analyse data of all monkeys together 

sel_area = {'PMd','PRR'};

sel_tuning = 1; % select units depending on tuning properties
% 1 - motor-goal tuned during late memory period
% 2 - mem-RT (same direction)
% 3 - mem-RT-move (same direction)

remove_Oc = 1;      % 1 - remove trials in which pre-cue occurs at unit's PD (to minimize visual activities)

interpolated = 0;   % 1 - interpolate activities | 0 - activities at 4 reach locations
normalized = 0;     % 1 - normalize activities to maximum in memory period | 0 - no normalization

% plot?
visualize = 1;
surf_plot = 0 & visualize;
trace_plot = 1 & visualize;
if surf_plot, interpolated = 1; normalized = 1; end % if surf plot, we need to interpolate and normalize
if trace_plot, normalized = 0; end

% set alignment modes - see presets in getSpd_hdf52mat.m
align_mode = {'preCue','go','move'};%,'touch' 'preCue','lateMem',

% set spd kernel (if other kernels are needed, go back and generate them in main_plx2hdf5.m)
spd_type = 'e50'; % e20 | e50 | g20 | g50

% call getSpd_hdf52mat to get data and save them in mat files
new_filename = getSpd_hdf52mat(data_folder,mat_folder,pool_monkey,sel_monkey,sel_area,align_mode,spd_type,sel_tuning);

%% go through the data files one by one
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load spd_table - analyses %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
relat_dir_name = {'PD','Op','OD','Om'}; % orthogonal cue/ noncue (CueDir)

for nf = 1:numel(new_filename)
    disp(['Loading spd data from: ',mat_folder new_filename{nf} '.mat']);
    load([mat_folder new_filename{nf} '.mat']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % some preprocessing: cue|bias direction shift & normalization %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % convert to categorical for easy handling
    spd_table.unitName = categorical(spd_table.unitName);
%     spd_table.condName = categorical(spd_table.condName);
    
    unit_prf_dir = unique(spd_table(:,1:4));
    n_unit_sel = height(unit_prf_dir);
    
    n_bin = length(time_bin);
    
    % get bias conditions
    varName_list = spd_table.Properties.VariableNames;
    if any(ismember(varName_list,'BiasType')), bias_name = 'BiasType'; end
    if any(ismember(varName_list,'BiasTypeAbs')), bias_name = 'BiasTypeAbs'; end
    
    dir_list = unique(spd_table.(setting.dtun));
    n_dir = length(dir_list);
    bias_list = unique(spd_table.(bias_name));
    n_bias = length(bias_list);
    cond_names = sorter_Varname;
    n_cond = length(cond_names);
    sorter_names = fieldnames(setting.sorter)';
    
    myMap = colormap('jet'); close(gcf);
    cols = {myMap(1,:), myMap(21,:), myMap(42,:), myMap(64,:)};
    
    
    if normalized
        suff = '_norm';
        fr_name = 'mean_spikeDensity_norm';
        yLim = 2;
    else
        suff = '';
        fr_name = 'mean_spikeDensity';
        yLim = 40;
    end
    
    if interpolated
        suff = [suff '_interp'];
    end
    
    dir_interp = (0:359)';
    n_dir_interp = 360;
    shift_pref_dir = 90;
    
    if remove_Oc
        % remove trials in which cueDir = maxDir
        remove_ind = spd_table.cueDir == spd_table.maxDir;
        spd_table = spd_table(~remove_ind,:);
        suff = [suff '_nonCue'];
    end
    % rotate bias direaction relative to PD - in variable: BiasDir_shift
    tmpDir_shift = spd_table.(setting.dtun);
    for uu = 1:n_unit_sel
        tmpAngle = circshift(dir_interp, shift_pref_dir - round(unit_prf_dir{uu,2}));
        cur_unit = spd_table.unitName == unit_prf_dir{uu,1};
        for dd = 1:n_dir
            cur_dir = spd_table.(setting.dtun) == dir_list(dd);
            tmpDir_shift(cur_unit&cur_dir) = find(tmpAngle == dir_list(dd));
        end
    end
    spd_table.([setting.dtun '_shift']) = (tmpDir_shift*2*pi)/360;
    
    tmp_relat_dir = spd_table.(setting.dtun);
    for dd = 1:n_dir
        cur_dir = spd_table.maxDir == dir_list(dd);
        rot_dir = circshift(dir_list,1-dd);
        for dc = 1:n_dir
            tmp_relat_dir(cur_dir & tmp_relat_dir == rot_dir(dc)) = dc;
        end
    end
    spd_table.([setting.dtun '_relat']) = relat_dir_name(tmp_relat_dir)';
    spd_table.([setting.dtun '_relat']) = categorical(spd_table.([setting.dtun '_relat']));
    
    % calculate normalized activities (max div) - in variable: spikeDensity_norm
    tmpspikeDensity_norm = spd_table.spikeDensity;
    for uu = 1:n_unit_sel
        cur_unit = spd_table.unitName == unit_prf_dir{uu,1};
        % normalize
        tmpspikeDensity_norm(cur_unit,:) = tmpspikeDensity_norm(cur_unit,:) ./ unique(spd_table.maxAct(cur_unit));
    end
    spd_table.spikeDensity_norm = tmpspikeDensity_norm;
    
    % add unit properties from spk analyses
    spk_path = strrep(result_path,'spd','spk');
    spk_filename = strrep(new_filename{nf},['spd' spd_type],'spk');
    spk_filename = strrep(spk_filename,num2str(sel_tuning),'');
    spk_filename = strrep(spk_filename,setting.name,'lateMem_bias');
    if exist([spk_path spk_filename '_graded_unit.mat'],'file')
        load([spk_path spk_filename '_graded_unit.mat'])
        tmp_graded_spk = spd_table.maxDir;
        tmp_graded = unit_prf_dir.maxDir;
        for uu = 1:n_unit_sel
            cur_unit = spd_table.unitName == unit_estimate.unitName(uu);
            sel_unit = unit_prf_dir.unitName == unit_estimate.unitName(uu);
            if unit_h.PD(uu) && unit_h.OD(uu), tmp_code = 3;
            elseif unit_h.PD(uu) && ~unit_h.OD(uu), tmp_code = 2;
            elseif ~unit_h.PD(uu) && unit_h.OD(uu), tmp_code = 1;
            else tmp_code = 0;
            end
            tmp_graded_spk(cur_unit) = tmp_code;
            tmp_graded(sel_unit) = tmp_code;
        end
        spd_table.gradedType = tmp_graded_spk;
        unit_prf_dir.gradedType = tmp_graded;
    else
        disp('no unit graded info - run spk_unit_graded')
    end
    
    
    % calculate mean activities for each unit - in table: spd_table_mean
    spd_mean = grpstats(spd_table,{'unitName',setting.dtun,'condName'},'mean','DataVars',{'spikeDensity','spikeDensity_norm'});
    spd_mean = unique(outerjoin(spd_mean,spd_table,'Keys',{'unitName',setting.dtun,'condName'},'Type','left','RightVariables',{'prefDir','maxAct','maxDir','gradedType'}),'rows');
    if interpolated
        if ~exist([result_path new_filename{nf} '_tuningProg' suff '.mat'],'file')
            % interpolate activities for tuning progression plots
            spd_tuning = nan(n_dir_interp,n_bin,n_cond,n_unit_sel);
            nUnit_cond = zeros(n_cond,1);
            unit_act = [];
            for uu = 1:n_unit_sel
                cur_unit = spd_mean.unitName == unit_prf_dir{uu,1};
                cur_PD = round(unit_prf_dir{uu,2});
                tmp_unit = spd_mean(cur_unit,:);
                for cc = 1:n_cond
                    cur_cond = tmp_unit.condName == cond_names(cc);
                    tmpComponents = nan(n_dir,n_bin);
                    for dd = 1:n_dir
                        cur_dir = tmp_unit.(setting.dtun) == dir_list(dd);
                        tmp_density = tmp_unit.(fr_name);
                        if ~isempty(tmp_density(cur_cond&cur_dir,:))
                            tmpComponents(dd,:) = tmp_density(cur_cond&cur_dir,:);
                        end
                    end
                    % interpolate
                    if ~any(isnan(tmpComponents))
                        nUnit_cond(cc) = nUnit_cond(cc) + 1;
                        
                        tmpTuning = interpft(tmpComponents,n_dir_interp,1);
                        tmpTuning = circshift(tmpTuning,shift_pref_dir - cur_PD);
                        spd_tuning(:,:,cc,uu) = tmpTuning;
                        
                        tmp_circ = circshift(dir_list,1-shift_prefDir);
                        for dd = 1:n_dir
                            tmp_ind = tmp_circ(dd)+1;
                            tmp_row = cell2table([{unit_prf_dir{uu,1}},cond_names(cc),relat_dir_name(dd),{tmpTuning(tmp_ind,:)}],...
                                'VariableNames',{'unitName','condName','gradedType',[setting.dtun '_relat'],fr_name});
                            unit_act = [unit_act;tmp_row];
                        end
                    end
                end
            end
            tuningProgressMat = imean(spd_tuning,4);
            
            cond_count = [cell2table(cond_names','VariableNames',{'condName'}),array2table(nUnit_cond,'VariableNames',{'GroupCount'})];
            cond_count.condName = categorical(cond_count.condName);
            cond_count.condName = reordercats(cond_count.condName,cond_names);
            cond_count = sortrows(cond_count,'condName');
            
            disp(['Saving tuningProg to: ',result_path new_filename{nf} '_tuningProg' suff '.mat']);
            save([result_path new_filename{nf} '_tuningProg' suff '.mat'],'unit_act','spd_tuning','cond_count', 'tuningProgressMat','setting', '-v7.3');
            
        else
            disp(['Loading tuningProg from: ',result_path new_filename{nf} '_tuningProg' suff '.mat']);
            load([result_path new_filename{nf} '_tuningProg' suff '.mat']);
        end
        
        if surfPlot
            temp_fig_path = checkDirectory(figure_path, 'tuningProg/',1);
            for cc = 1:n_cond
                h = figure;
                surf(tuningProgressMat(:,:,cc)');
                grid on
                shading interp
                colormap jet
                view(62, 54)
                
                title({[new_filename{nf} suff]; ['(n=' num2str(cond_count.GroupCount(cc)) ')']},'interpreter', 'none');
                xlabel('angle (deg) rel. to pref. direction','FontWeight','bold');
                set(gca,'XTick',1:90:361);
                
                % mark PD in graph (for 0 and 90 degrees only)
                if shift_prefDir == 0
                    set(gca,'XTickLabel',{'0 (PD)','90','180','270','360'});
                elseif shift_prefDir == 90
                    set(gca,'XTickLabel',{'0','90 (PD)','180','270','360'});
                elseif shift_prefDir == 180
                    set(gca,'XTickLabel',{'0','90','180 (PD)','270','360'});
                elseif shift_prefDir == 270
                    set(gca,'XTickLabel',{'0','90','180','270 (PD)','360'});
                end
                xlim([0 n_dir_interp]);
                ylabel(['time (ms) rel. to: ' align{1} ],'FontWeight','bold');
                set(gca,'YTick',1:100:n_bin);
                set(gca,'YTickLabel',time_bin(1):100:time_bin(end));
                ylim([0 n_bin]);
                zlim([0 yLim]);
                %             set(gca,'CLim',[0.5 1])
                colorbar
                saveas(h,[temp_fig_path new_filename{nf} '_tuningProg_' cond_names{cc} suff '.fig']);
                save2pdf([temp_fig_path new_filename{nf} '_tuningProg_' cond_names{cc} suff '.pdf'],h);
                caf
            end
        end
        
    else
        if trace_plot
            if ~exist([result_path new_filename{nf} '_unitAct' suff '.mat'],'file')
                unit_act = spd_mean(:,{'unitName','condName','gradedType',fr_name});
                tmp_relat_dir = spd_mean.(setting.dtun);
                for dd = 1:n_dir
                    cur_dir = spd_mean.maxDir == dir_list(dd);
                    rot_dir = circshift(dir_list,1-dd);
                    for dc = 1:n_dir
                        tmp_relat_dir(cur_dir & tmp_relat_dir == rot_dir(dc)) = dc;
                    end
                end
                unit_act.([setting.dtun '_relat']) = relat_dir_name(tmp_relat_dir)';
                
                unit_unstack = unstack(unit_act,fr_name,[setting.dtun '_relat']);
                unit_mean = grpstats(unit_unstack,{'condName'},'mean','DataVars',{'PD'});
                unit_mean.Properties.RowNames = {};
                cond_count = unit_mean(:,{'condName','GroupCount'});
                cond_count.Properties.RowNames = {};
                if height(cond_count) < n_cond
                    ind_missing = ~ismember(cond_names,cond_count.condName);
                    tmp_condName = cell2table([cellstr(cond_count.condName);cond_names(ind_missing)],'VariableNames',{'condName'});
                    tmp_GroupCount = array2table([cond_count.GroupCount;zeros(sum(ind_missing),1)],'VariableNames',{'GroupCount'});
                    cond_count = [tmp_condName,tmp_GroupCount];
                    cond_count.condName = categorical(cond_count.condName);
                    cond_count.condName = reordercats(cond_count.condName,cond_names);
                    cond_count = sortrows(cond_count,'condName');
                end
                disp(['Saving unit activity to: ',result_path new_filename{nf} '_unitAct' suff '.mat']);
                save([result_path new_filename{nf} '_unitAct' suff '.mat'],'unit_act','cond_count','setting');
            else
                disp(['Loading unit activity from: ',result_path new_filename{nf} '_unitAct' suff '.mat']);
                load([result_path new_filename{nf} '_unitAct' suff '.mat']);
            end
        end
    end
    if trace_plot
        temp_fig_path = checkDirectory(figure_path, 'tracePlots/',1);
        sorter_list = fieldnames(setting.sorter);
        sub_sorter_list = cell(length(sorter_list),1);
        
        unit_trace = unit_act;
        tmp_condName_split = cell2table(table2array(cell2table(cellfun(@(x) strsplit(x,'_'),cellstr(unit_trace.condName),'UniformOutput',0))),'VariableNames',sorter_list);
        unit_trace = [unit_trace,tmp_condName_split];
        unit_trace.([setting.dtun '_relat']) = categorical(unit_trace.([setting.dtun '_relat']));
        
        ind_bias = ismember(sorter_list,bias_name);
        n_bias = length(unique(tmp_condName_split{:,ind_bias}));
        
        tmp_condName_split = cell2table(table2array(cell2table(cellfun(@(x) strsplit(x,'_'),cellstr(cond_count.condName),'UniformOutput',0))),'VariableNames',sorter_list);
        if any(~ind_bias)
            cond_trace = unique(tmp_condName_split(:,~ind_bias),'rows');
            cond_list = cond_trace.Properties.VariableNames;
            unit_count = zeros(height(cond_trace),n_bias);
            for cc = 1:height(cond_trace)
                tmp_ind = zeros(height(cond_count),length(cond_list));
                tmp_sel = zeros(height(unit_trace),length(cond_list));
                fig_suff = '';
                for ss = 1:length(cond_list)
                    tmp_ind(:,ss) = ismember(tmp_condName_split.(cond_list{ss}),cond_trace.(cond_list{ss})(cc));
                    tmp_sel(:,ss) = ismember(unit_trace.(cond_list{ss}),cond_trace.(cond_list{ss})(cc));
                    fig_suff = [fig_suff '_' cond_trace.(cond_list{ss}){cc}];
                end
                cur_cond = all(tmp_ind,2);
                unit_count(cc,:) = cond_count.GroupCount(cur_cond);
                sel_dir = (unit_trace.([setting.dtun '_relat']) == 'PD' | unit_trace.([setting.dtun '_relat']) == 'OD');
                cur_sel = all(tmp_sel,2) & sel_dir;
                
                g = gramm('x',repmat(time_bin,height(unit_trace),1),'y',unit_trace.(fr_name),...
                    'lightness',unit_trace.(bias_name),'linestyle',unit_trace.([setting.dtun '_relat']),...
                    'subset',cur_sel,'color',unit_trace.([setting.dtun '_relat']));
                g.stat_summary('Type','bootci')
                g.set_order_options('color',{'PD','OD'},'linestyle',{'PD','OD'})
                g.set_names('x',['time rel. to' align],'y', fr_name)
                g.axe_property('ylim',[0 yLim])
                g.geom_vline('xintercept',0,'style','k:')
                g.set_title([new_filename{nf} fig_suff suff ' (n = ' num2str(unit_count(cc,:)) ')'])
                h = figure('Name',[new_filename{nf} fig_suff suff],'Position',get(0,'ScreenSize'));
                g.draw()
                saveas(h,[temp_fig_path new_filename{nf} '_trace_' fig_suff suff '.fig']);
                save2pdf([temp_fig_path new_filename{nf} '_trace_' fig_suff suff '.pdf'],h);
                caf
                
                g = gramm('x',repmat(time_bin,height(unit_trace),1),'y',unit_trace.(fr_name),...
                    'lightness',unit_trace.(bias_name),'linestyle',unit_trace.([setting.dtun '_relat']),...
                    'subset',all(tmp_sel,2),'color',unit_trace.([setting.dtun '_relat']));
                g.stat_summary('Type','bootci')
                g.set_order_options('color',{'PD','OD','Op','Om'},'linestyle',{'PD','OD','Op','Om'})
                g.set_names('x',['time rel. to' align],'y', fr_name)
                g.axe_property('ylim',[0 yLim])
                g.geom_vline('xintercept',0,'style','k:')
                g.set_title([new_filename{nf} fig_suff suff ' (n = ' num2str(unit_count(cc,:)) ')'])
                h = figure('Name',[new_filename{nf} fig_suff suff],'Position',get(0,'ScreenSize'));
                g.draw()
                saveas(h,[temp_fig_path new_filename{nf} '_trace_allDir' fig_suff suff '.fig']);
                save2pdf([temp_fig_path new_filename{nf} '_trace_allDir' fig_suff suff '.pdf'],h);
                caf
                
                gradedType_list = unique(unit_trace.gradedType);
                unit_trace.(bias_name) = categorical(unit_trace.(bias_name));
                for gg = 1:length(gradedType_list)
                    graded_sel = unit_trace.gradedType == gradedType_list(gg);
                    n_sel = zeros(1,n_bias);
                    fig_suff_g = [fig_suff '_selGraded' num2str(gradedType_list(gg))];
                    for bb = 1:n_bias
                        n_sel(bb) = sum(cur_sel&graded_sel&unit_trace.(bias_name)==[bias_name num2str(bias_list(bb))]);
                    end
                    g = gramm('x',repmat(time_bin,height(unit_trace),1),'y',unit_trace.(fr_name),...
                        'lightness',cellstr(unit_trace.(bias_name)),'linestyle',unit_trace.([setting.dtun '_relat']),...
                        'subset',cur_sel&graded_sel,'color',unit_trace.([setting.dtun '_relat']));
                    g.stat_summary('Type','bootci')
                    g.set_order_options('color',{'PD','OD'},'linestyle',{'PD','OD'})
                    g.set_names('x',['time rel. to' align],'y', fr_name)
                    g.axe_property('ylim',[0 yLim])
                    g.geom_vline('xintercept',0,'style','k:')
                    g.set_title([new_filename{nf} fig_suff_g suff ' (n = ' num2str(n_sel) ')'])
                    h = figure('Name',[new_filename{nf} fig_suff_g suff],'Position',get(0,'ScreenSize'));
                    g.draw()
                    saveas(h,[temp_fig_path new_filename{nf} '_trace_' fig_suff_g suff '.fig']);
                    save2pdf([temp_fig_path new_filename{nf} '_trace_' fig_suff_g suff '.pdf'],h);
                    caf
                end
            end
        else
            cur_sel = (unit_trace.([setting.dtun '_relat']) == 'PD' | unit_trace.([setting.dtun '_relat']) == 'OD');
            g = gramm('x',repmat(time_bin,height(unit_trace),1),'y',unit_trace.(fr_name),...
                'lightness',unit_trace.(bias_name),'linestyle',unit_trace.([setting.dtun '_relat']),...
                'subset',cur_sel,'color',unit_trace.([setting.dtun '_relat']));
            g.stat_summary('Type','bootci')
            g.set_order_options('color',{'PD','OD'},'linestyle',{'PD','OD'})
            g.set_names('x',['time rel. to' align],'y', fr_name)
            g.axe_property('ylim',[0 yLim])
            g.geom_vline('xintercept',0,'style','k:')
            g.set_title([new_filename{nf} suff ' (n = ' num2str(cond_count.GroupCount') ')'])
            h = figure('Name',[new_filename{nf} suff],'Position',get(0,'ScreenSize'));
            g.draw()
            saveas(h,[temp_fig_path new_filename{nf} '_trace_' suff '.fig']);
            save2pdf([temp_fig_path new_filename{nf} '_trace_' suff '.pdf'],h);
            caf
            
            g = gramm('x',repmat(time_bin,height(unit_trace),1),'y',unit_trace.(fr_name),...
                'lightness',unit_trace.(bias_name),'linestyle',unit_trace.([setting.dtun '_relat']),...
                'color',unit_trace.([setting.dtun '_relat'])); %'subset',cur_sel
            g.stat_summary('Type','bootci')
            g.set_order_options('color',{'PD','OD','Op','Om'},'linestyle',{'PD','OD','Op','Om'})
            g.set_names('x',['time rel. to' align],'y', fr_name)
            g.axe_property('ylim',[0 yLim])
            g.geom_vline('xintercept',0,'style','k:')
            g.set_title([new_filename{nf} suff ' (n = ' num2str(cond_count.GroupCount') ')'])
            h = figure('Name',[new_filename{nf} suff],'Position',get(0,'ScreenSize'));
            g.draw()
            saveas(h,[temp_fig_path new_filename{nf} '_trace_allDir' suff '.fig']);
            save2pdf([temp_fig_path new_filename{nf} '_trace_allDir' suff '.pdf'],h);
            caf
            
            gradedType_list = unique(unit_trace.gradedType);
            unit_trace.(bias_name) = categorical(unit_trace.(bias_name));
            for gg = 1:length(gradedType_list)
                graded_sel = unit_trace.gradedType == gradedType_list(gg);
                n_sel = zeros(1,n_bias);
                fig_suff = ['_selGraded' num2str(gradedType_list(gg))];
                for bb = 1:n_bias
                    n_sel(bb) = sum(cur_sel&graded_sel&unit_trace.(bias_name)==[bias_name num2str(bias_list(bb))]);
                end
                g = gramm('x',repmat(time_bin,height(unit_trace),1),'y',unit_trace.(fr_name),...
                    'lightness',cellstr(unit_trace.(bias_name)),'linestyle',unit_trace.([setting.dtun '_relat']),...
                    'subset',cur_sel&graded_sel,'color',unit_trace.([setting.dtun '_relat']));
                g.stat_summary('Type','bootci')
                g.set_order_options('color',{'PD','OD'},'linestyle',{'PD','OD'})
                g.set_names('x',['time rel. to' align],'y', fr_name)
                g.axe_property('ylim',[0 yLim])
                g.geom_vline('xintercept',0,'style','k:')
                g.set_title([new_filename{nf} '_trace' fig_suff suff ' (n = ' num2str(n_sel) ')'])
                h = figure('Name',[new_filename{nf} '_trace' fig_suff suff],'Position',get(0,'ScreenSize'));
                g.draw()
                saveas(h,[temp_fig_path new_filename{nf} '_trace' fig_suff suff '.fig']);
                save2pdf([temp_fig_path new_filename{nf} '_trace' fig_suff suff '.pdf'],h);
                caf
            end
        end
        
        % time diff plot
        unit_trace_unstack = unstack(unit_trace,fr_name,([setting.dtun '_relat']));
        unit_trace_unstack.diff = unit_trace_unstack.PD - unit_trace_unstack.OD;
        
        if any(~ind_bias)
            for cc = 1:height(cond_trace)
                tmp_sel = zeros(height(unit_trace_unstack),length(cond_list));
                fig_suff = '';
                for ss = 1:length(cond_list)
                    tmp_sel(:,ss) = ismember(unit_trace_unstack.(cond_list{ss}),cond_trace.(cond_list{ss})(cc));
                    fig_suff = [fig_suff '_' cond_trace.(cond_list{ss}){cc}];
                end
                cur_sel = all(tmp_sel,2);
                
                g = gramm('x',repmat(time_bin,height(unit_trace_unstack),1),'y',unit_trace_unstack.diff,'lightness',unit_trace_unstack.(bias_name),'subset',cur_sel);
                g.stat_summary('Type','bootci')
                g.set_names('x',['time rel. to' align],'y', [fr_name '_diffOD'])
                g.axe_property('ylim',[-yLim*3/4 yLim*3/4])
                g.geom_vline('xintercept',0,'style','k:')
                g.geom_hline('yintercept',0,'style','k:')
                g.set_title([new_filename{nf} '_diffOD' fig_suff suff ' (n = ' num2str(unit_count(cc,:)) ')'])
                h = figure('Name',[new_filename{nf} '_diffOD' fig_suff suff],'Position',get(0,'ScreenSize'));
                g.draw()
                saveas(h,[temp_fig_path new_filename{nf} '_diffTraceOD_' fig_suff suff '.fig']);
                save2pdf([temp_fig_path new_filename{nf} '_diffTraceOD_' fig_suff suff '.pdf'],h);
                caf
            end
        else
            g = gramm('x',repmat(time_bin,height(unit_trace_unstack),1),'y',unit_trace_unstack.diff,'lightness',unit_trace_unstack.(bias_name));
            g.stat_summary('Type','bootci')
            g.set_names('x',['time rel. to' align],'y', [fr_name '_diffOD'])
            g.axe_property('ylim',[-yLim*3/4 yLim*3/4])
            g.geom_vline('xintercept',0,'style','k:')
            g.geom_hline('yintercept',0,'style','k:')
            g.set_title([new_filename{nf} suff ' (n = ' num2str(cond_count.GroupCount') ')'])
            h = figure('Name',[new_filename{nf} '_diffOD' suff],'Position',get(0,'ScreenSize'));
            g.draw()
            saveas(h,[temp_fig_path new_filename{nf} '_diffTraceOD_' suff '.fig']);
            save2pdf([temp_fig_path new_filename{nf} '_diffTraceOD_' suff '.pdf'],h);
            caf
        end
        
        % time diff plot2 (diff from orth)
        for bn = 1:n_bin
            unit_trace_unstack.orth(:,bn) = nanmean([unit_trace_unstack.Om(:,bn),unit_trace_unstack.Op(:,bn)],2);
        end
        unit_trace_unstack.diff = unit_trace_unstack.PD - unit_trace_unstack.orth;
        
        if any(~ind_bias)
            for cc = 1:height(cond_trace)
                tmp_sel = zeros(height(unit_trace_unstack),length(cond_list));
                fig_suff = '';
                for ss = 1:length(cond_list)
                    tmp_sel(:,ss) = ismember(unit_trace_unstack.(cond_list{ss}),cond_trace.(cond_list{ss})(cc));
                    fig_suff = [fig_suff '_' cond_trace.(cond_list{ss}){cc}];
                end
                cur_sel = all(tmp_sel,2);
                
                g = gramm('x',repmat(time_bin,height(unit_trace_unstack),1),'y',unit_trace_unstack.diff,'lightness',unit_trace_unstack.(bias_name),'subset',cur_sel);
                g.stat_summary('Type','bootci')
                g.set_names('x',['time rel. to' align],'y', [fr_name '_diffOrth'])
                g.axe_property('ylim',[-yLim*3/4 yLim*3/4])
                g.geom_vline('xintercept',0,'style','k:')
                g.geom_hline('yintercept',0,'style','k:')
                g.set_title([new_filename{nf} '_diffOrth' fig_suff suff ' (n = ' num2str(unit_count(cc,:)) ')'])
                h = figure('Name',[new_filename{nf} '_diffOrth' fig_suff suff],'Position',get(0,'ScreenSize'));
                g.draw()
                saveas(h,[temp_fig_path new_filename{nf} '_diffTraceOrth_' fig_suff suff '.fig']);
                save2pdf([temp_fig_path new_filename{nf} '_diffTraceOrth_' fig_suff suff '.pdf'],h);
                caf
            end
        else
            g = gramm('x',repmat(time_bin,height(unit_trace_unstack),1),'y',unit_trace_unstack.diff,'lightness',unit_trace_unstack.(bias_name));
            g.stat_summary('Type','bootci')
            g.set_names('x',['time rel. to' align],'y', [fr_name '_diffOrth'])
            g.axe_property('ylim',[-yLim*3/4 yLim*3/4])
            g.geom_vline('xintercept',0,'style','k:')
            g.geom_hline('yintercept',0,'style','k:')
            g.set_title([new_filename{nf} suff ' (n = ' num2str(cond_count.GroupCount') ')'])
            h = figure('Name',[new_filename{nf} '_diffOrth' suff],'Position',get(0,'ScreenSize'));
            g.draw()
            saveas(h,[temp_fig_path new_filename{nf} '_diffTraceOrth_' suff '.fig']);
            save2pdf([temp_fig_path new_filename{nf} '_diffTraceOrth_' suff '.pdf'],h);
            caf
        end
    end
    if graded_unit_latency
        varName_list = spd_table.Properties.VariableNames;
        if any(ismember(varName_list,'BiasType')), biasName = 'BiasType'; end
        if any(ismember(varName_list,'BiasTypeAbs')), biasName = 'BiasTypeAbs'; end
        
        Bias_ratio = [0.5 0.375 0.25 0.125; 0.5 0.625 0.75 0.875];
        nB = size(Bias_ratio,2);
        bias_strength = zeros(1,nB);
        bias_range = zeros(1,nB*2-1);
        for b = 1:nB
            bias_strength(b) = (Bias_ratio(2,b)-Bias_ratio(1,b))/(Bias_ratio(2,b)+Bias_ratio(1,b));
            bias_range(b+3) = (Bias_ratio(2,b)-Bias_ratio(1,b))/(Bias_ratio(2,b)+Bias_ratio(1,b));
            bias_range(5-b) = (Bias_ratio(1,b)-Bias_ratio(2,b))/(Bias_ratio(2,b)+Bias_ratio(1,b));
        end
        switch biasName
            case 'BiasType',    spd_table.BiasLevel = bias_range(spd_table.(biasName)+1)';
            case 'BiasTypeAbs', spd_table.BiasLevel = bias_strength(spd_table.(biasName)+1)';
        end
        temp_fig_path = checkDirectory(figure_path, 'graded_unit_latency/',1);
        relat_dir_name2 = {'PD','OD'};
        if ~exist([result_path new_filename{nf} '_graded_unit_latency.mat'],'file')
            graded_unit = [];
            gradedArray_unit = [];
            tmp_PDOD = spd_table.([setting.dtun '_relat']) == 'PD' | spd_table.([setting.dtun '_relat']) == 'OD';
            for uu = 1:nUnit_sel
                tmp_u = unit_prf_dir{uu,1};
                tmp_unit = spd_table(spd_table.unitName == tmp_u & tmp_PDOD,:);
                for dd = 1:nDir/2
                    % rearrange categories to get 'PD', 'OD', etc. as the first
                    % element when sorted
                    tmp_dir_name = circshift(relat_dir_name2',1-dd);
                    tmp_unit.([setting.dtun '_relat']) = categorical(tmp_unit.([setting.dtun '_relat']),...
                        tmp_dir_name,'Ordinal',true);
                    tmp_unit = sortrows(tmp_unit,[setting.dtun '_relat']);
                    tmp_unit.([setting.dtun '_relat']) = categorical(tmp_unit.([setting.dtun '_relat']),...
                        tmp_dir_name,'Ordinal',false);
                    
                    
                    estimate = nan(1,n_bin);
                    pValue = nan(1,n_bin);
                    for bb = 1:n_bin
                        % glm fit
                        tmp_unit.curr_bin = tmp_unit.spikeDensity_norm(:,bb);
                        glm_graded_unit = fitglm(tmp_unit,['curr_bin ~ BiasLevel * ' [setting.dtun '_relat'] ' ']);
                        tmp_ci = coefCI(glm_graded_unit);
                        tmp_ind = strcmp(glm_graded_unit.CoefficientNames,'BiasLevel');
                        estimate(bb) = glm_graded_unit.Coefficients{tmp_ind,'Estimate'};
                        pValue(bb) = glm_graded_unit.Coefficients{tmp_ind,'pValue'};
                    end
                    tmp_graded = [unit_prf_dir(uu,1),table(estimate),table(pValue),table(tmp_dir_name(1),'VariableNames',{'relatDir'})];
                        graded_unit = [graded_unit;tmp_graded];
                end
            end
            graded_unit.relatDir = categorical(graded_unit.relatDir);
            graded_unit.h = graded_unit.pValue < 0.05;
            
            unit_estimate = unstack(graded_unit(:,{'unitName','estimate','relatDir'}),'Estimate','relatDir');
            unit_h = unstack(graded_unit(:,{'unitName','h','relatDir'}),'h','relatDir');
            
            graded_name = {'gradedBoth','gradedPrf','gradedOpp','notGraded'};
            graded_unit_sum = cell(size(unit_h.PD));
            graded_unit_sum(unit_h.PD & unit_h.OD) = graded_name(1);
            graded_unit_sum(unit_h.PD & ~unit_h.OD) = graded_name(2);
            graded_unit_sum(~unit_h.PD & unit_h.OD) = graded_name(3);
            graded_unit_sum(~unit_h.PD & ~unit_h.OD) = graded_name(4);
            graded_unit_num = sum([unit_h.PD & unit_h.OD,unit_h.PD & ~unit_h.OD,~unit_h.PD & unit_h.OD,~unit_h.PD & ~unit_h.OD]);
            
            save([result_path new_filename{nf} '_graded_unit_latency.mat'],'unit_estimate','graded_unit_num','graded_unit_sum','unit_h')
        else
            load([result_path new_filename{nf} '_graded_unit_latency.mat'])
        end
        
    end
end

