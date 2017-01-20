%%%%%%%%%%%%%
% set paths %
%%%%%%%%%%%%%
if ismac
    data_folder = '/Volumes/Passport/dataAnalysis/data/';
    result_path = '/dataAnalysis/results/spk/';
else % in case used on Windows
    result_path = 'D:\results\hdf5\spk\';
    data_folder = 'D:\data\Decision_bias\';
end
mat_folder = '/dataAnalysis/mat/spk/';
% mat_folder = [data_folder 'mat' filesep 'spd' filesep];
figure_path = [result_path 'figures' filesep];

if ~exist(figure_path,'dir'), mkdir(figure_path); end
if ~exist(mat_folder,'dir'), mkdir(mat_folder); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select monkeys, areas, tasks %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sel_area = {'PMd','PRR'};
area_comp = 0 & numel(sel_area)>1;  

sel_monkey = {'Hum','Kas'};
pool_monkey = 1 & numel(sel_monkey)>1;  
% 0: analyse data of all monkeys separately
% 1: analyse data of all monkeys together 

sel_taskVar = {'rew_rnd','rew_HL'};
pool_taskVar = 1 & numel(sel_taskVar)>1;
% 0: analyse taskVar separately
% 1: analyse taskVar together 

remove_Oc = 0;
normInterp_max = 0;

graded_pop = 0;
test_graded_unit = 1;

show_unsmooth = 0;
if show_unsmooth, normInterp_max = 1; end
vonMises_fit = 0;
if vonMises_fit, remove_Oc = 1; graded_pop = 1; end

card_oblique = 0;

glmm_test = 0;
parallel_sig = 0;
parallel_noise = 0;
single_corr_plot = 0;
cardOnly = 0;
corr_selDir = 0;

test_ROC = 0;
single_ROC_plot = 0;
decoding = 0;

%%%%%%%%%%%%%%%%%%%%%%%
% set alignment modes %
%%%%%%%%%%%%%%%%%%%%%%%

align_mode = {'lateMem_bias'};

new_filename = getSpk_hdf52mat(data_folder,mat_folder,sel_monkey,pool_monkey,sel_area,sel_taskVar,pool_taskVar,align_mode);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load spk_table - analyses %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nf = 1:numel(new_filename)
    disp(['Loading spk data from: ',mat_folder new_filename{nf} '.mat']);
    load([mat_folder new_filename{nf} '.mat']);
    
    % clean-up
    spk_table(spk_table.moveDir == spk_table.cueDir,:) = [];
    
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % some preprocessing: cue|bias direction shift & normalization %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % convert unitName to categorical for easy handling
    spk_table.unitName = categorical(spk_table.unitName);
    
    % get bias conditions
    varName_list = spk_table.Properties.VariableNames;
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
    

    % get PrefDir
    unitPrfDir = unique(spk_table(:,1:4));
    nUnit_sel = height(unitPrfDir);
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % neural tuning codes %
    %%%%%%%%%%%%%%%%%%%%%%%
    % get unit summary (run main_tuningSummary)
    if pool_monkey, cur_monkey = 'bothMonkeys'; end
    if strfind(new_filename{nf},'PMd'), cur_area = 'PMd'; 
    elseif strfind(new_filename{nf},'PRR'), cur_area = 'PRR'; 
    end
    summary_path = ['/dataAnalysis/results/tuningSummary/unit_categories_' cur_area '_' cur_monkey '.mat'];
    if exist(summary_path,'file')
        load(summary_path);
        tmp_tune = spk_table.maxDir;
        unit_tune = unitPrfDir.maxDir;
        for uu = 1:nUnit_sel
            cur_unit_spk = spk_table.unitName == unitPrfDir{uu,1};
            cur_unit_tun = cur_tuning.unit_name == unitPrfDir{uu,1};
            tmp_tune(cur_unit_spk) = cur_tuning.tuneCode(cur_unit_tun);
            unit_tune(uu) = cur_tuning.tuneCode(cur_unit_tun);
        end
        spk_table.tune_code = tmp_tune;
        unitPrfDir.tune_code = unit_tune;
    end
    
    % get some condition & direction lists
    dir_list = unique(spk_table.(setting.dtun));
    nDir = length(dir_list);
    cond_list = unique(spk_table.(biasName));
    cond_name = unique(spk_table.condName);
    ncond = length(cond_list);
    
    % angle for activity interpolation
    dirInterp = (0:359)';
    nDirInterp = 360;
    if strfind(setting.dtun,'Cue')
        shift_dir = 180;
        x_label = {'-pi','-pi/2','0 (PD)','pi/2','pi'};
    else
        shift_dir = 90;
        x_label = {'-pi/2','0 (PD)','pi/2','pi','3pi/2'};
    end
    shift_rad = shift_dir * pi /180;
    
    suff = '';
    if remove_Oc
        suff = '_nonCue';
        % remove trials in which cueDir = maxDir
        remove_ind = spk_table.cueDir == spk_table.maxDir;
        spk_table = spk_table(~remove_ind,:);
        
        % when remove cued-orth direction, only one Orth (Op) in zero bias cond
        for du = 1:nDir
            rot_dir = circshift(dir_list,1-du)'; % rot_dir(1) = maxDir
            % select all units of the same maxDir
            ind_selUnit = spk_table.maxDir == dir_list(du);
            % in case balance cue appears on the other side
            ind_balOp = spk_table.BiasTypeAbs == 0 & spk_table.cueDir == rot_dir(3); % rot_dir(3) = opposite direction from maxDir
            ind_biasZeroOp = find(ind_selUnit & ind_balOp);
            n_biasZeroOp = sum(ind_selUnit & ind_balOp);
            ind_rand = randperm(n_biasZeroOp);
            ind_randHalf = ind_rand(1:round(n_biasZeroOp/2));
            % set bias direction to the other direction in half of the trials
            spk_table.(setting.dtun)(ind_biasZeroOp(ind_randHalf)) = repmat(rot_dir(4),length(ind_randHalf),1);
        end
    end
    
    % rotate angles to be relative to preferred direction
    tmpDir_shift = spk_table.(setting.dtun);
    for uu = 1:nUnit_sel
        cur_unit = spk_table.unitName == unitPrfDir{uu,1};
        tmpAngle = circshift(dirInterp, shift_dir - round(unitPrfDir{uu,2}) - 1);
        for dd = 1:nDir
            cur_dir = spk_table.(setting.dtun) == dir_list(dd);
            tmpDir_shift(cur_unit&cur_dir) = find(tmpAngle == dir_list(dd));
        end
    end
    spk_table.([setting.dtun '_shift']) = (tmpDir_shift*2*pi)/360; 

    % categorize cue|bias directions into PD (maxDir), +90, OD & -90
    relat_dir_name = {'PD','Op','OD','Om'};
    tmp_spk_table = [];
    for du = 1:nDir
        ind_selUnit = spk_table.maxDir == dir_list(du);
        tmp_sel_unit = spk_table(ind_selUnit,:);
        rot_dir = circshift(dir_list,1-du)';
        tmp_dir_relat = tmp_sel_unit.(setting.dtun);
        for dc = 1:nDir
            tmp_dir_relat(tmp_dir_relat == rot_dir(dc)) = dc;
        end
        tmp_sel_unit.([setting.dtun '_relat']) = tmp_dir_relat;
        tmp_spk_table = [tmp_spk_table;tmp_sel_unit];
    end
    tmp_spk_table.([setting.dtun '_relat_name']) = categorical(relat_dir_name(tmp_spk_table.([setting.dtun '_relat']))');
    spk_table = tmp_spk_table;
    
    if normInterp_max
        suff = [suff '_normInterpMax'];
        if ~exist([mat_folder new_filename{nf} suff '.mat'],'file')
            if strfind(setting.name,'lateMem_bias')
                spk_mean = grpstats(spk_table,{'unitName',biasName,setting.dtun,[setting.dtun '_shift']},'mean','DataVars',{'firingRate'});
                spk_mean = unique(outerjoin(spk_mean,spk_table,'Keys',{'unitName',biasName,setting.dtun},'Type','left','RightVariables',{'prefDir','maxAct','maxDir','condName'}),'rows');
                % take only extreme bias
                sel_cond = spk_mean.(biasName) == 3;
                unitInterpMax = unitPrfDir;
                unitInterpMax.Properties.VariableNames{'prefDir'} = 'interpMaxAct';
                for uu = 1:nUnit_sel
                    cur_unit = spk_mean.unitName == unitPrfDir{uu,1};
                    cur_PD = round(unitPrfDir{uu,2});
                    
                    tmp_unit = spk_mean(cur_unit&sel_cond,{biasName,setting.dtun,'mean_firingRate'});
                    tmpComponents = tmp_unit.mean_firingRate';
                    tmpTuning = interpft(tmpComponents,nDirInterp,2);
                    tmpTuning = circshift(tmpTuning,shift_dir - cur_PD,2);
                    unitInterpMax.interpMaxAct(uu) = tmpTuning(shift_dir);
                end
                save([mat_folder new_filename{nf} suff '.mat'], 'unitInterpMax')
            else
                disp('setting does not allow to compute interpolated maximum! - choose setting.name = lateMem_bias')
            end
        else
            load([mat_folder new_filename{nf} suff '.mat'])
        end
        
        % add to spk_table
        spk_table.interpMaxAct = spk_table.maxAct;
        for uu = 1:nUnit_sel
            cur_unit_spk_table = spk_table.unitName == unitInterpMax{uu,1};
            spk_table.interpMaxAct(cur_unit_spk_table) = unitInterpMax{uu,2};
        end
        
        % normalize to interpolated maximal activity
        spk_table.firingRate_norm = spk_table.firingRate./spk_table.interpMaxAct;
    else
        % normalize to maximal activity
        spk_table.firingRate_norm = spk_table.firingRate./spk_table.maxAct;
    end
  
       
    switch biasName
        case 'BiasType',    spk_table.BiasLevel = bias_range(spk_table.(biasName)+1)';
        case 'BiasTypeAbs', spk_table.BiasLevel = bias_strength(spk_table.(biasName)+1)';
    end
    
    
    % mark cardinal neurons : 45 deg around cardinal directions
    tmp_card = zeros(size(spk_table.prefDir));
    for du = 1:nDir
        if du == 1
            tmp_card = (abs(spk_table.prefDir - 0) < 22.5) | (abs(spk_table.prefDir - 360) < 22.5);
        else
            tmp_card = tmp_card | abs(spk_table.prefDir - dir_list(du)) < 22.5;
        end
    end
    spk_table.card = tmp_card;
    n_card = length(unique(spk_table.unitName(spk_table.card == 1)));
    
    unitPrfDir.card = ismember(unitPrfDir.unitName,unique(spk_table.unitName(spk_table.card == 1)));
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calculate & visualize mean FRs %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    spk_mean = grpstats(spk_table,{'unitName',biasName,setting.dtun,[setting.dtun '_shift']},'mean','DataVars',{'firingRate','firingRate_norm'});
    spk_mean = unique(outerjoin(spk_mean,spk_table,'Keys',{'unitName',biasName,setting.dtun},'Type','left','RightVariables',{'prefDir','maxAct','maxDir','condName'}),'rows');
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GRADED TEST: single units %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    temp_fig_path = checkDirectory(figure_path, 'graded_unit/',1);
    relat_dir_name2 = {'PD','OD'};
    if ~exist([result_path new_filename{nf} '_graded_unit.mat'],'file')
        graded_unit = [];
        gradedArray_unit = [];
        for uu = 1:nUnit_sel
            tmp_u = unitPrfDir{uu,1};
            tmp_unit = spk_table(spk_table.unitName == tmp_u,:);
            for dd = 1:nDir
                % rearrange categories to get 'PD', 'OD', etc. as the first
                % element when sorted
                tmp_dir_name = circshift(relat_dir_name',1-dd);
                tmp_unit.([setting.dtun '_relat_name']) = categorical(tmp_unit.([setting.dtun '_relat_name']),...
                    tmp_dir_name,'Ordinal',true);
                tmp_unit = sortrows(tmp_unit,[setting.dtun '_relat_name']);
                tmp_unit.([setting.dtun '_relat_name']) = categorical(tmp_unit.([setting.dtun '_relat_name']),...
                    tmp_dir_name,'Ordinal',false);
                % glm fit
                glm_graded_unit = fitglm(tmp_unit,['firingRate_norm ~ BiasLevel * ' [setting.dtun '_relat_name'] ' ']);
                tmp_ci = coefCI(glm_graded_unit);
                tmp_ind = strcmp(glm_graded_unit.CoefficientNames,'BiasLevel');
                tmp_graded = [unitPrfDir(uu,1),...
                    cell2table(tmp_dir_name(1),'VariableNames',{'relatDir'}),...
                    glm_graded_unit.Coefficients(tmp_ind,:),...
                    array2table(tmp_ci(tmp_ind,:),'VariableNames',{'lw_CI','up_CI'})];
                tmp_graded.Properties.RowNames = {};
                graded_unit = [graded_unit;tmp_graded];
            end    
        end
        graded_unit.relatDir = categorical(graded_unit.relatDir);
        graded_unit.h = graded_unit.pValue < 0.05;
    
        unit_estimate = unstack(graded_unit(:,{'unitName','Estimate','relatDir'}),'Estimate','relatDir');
        unit_h = unstack(graded_unit(:,{'unitName','h','relatDir'}),'h','relatDir');
    
        graded_name = {'gradedBoth','gradedPrf','gradedOpp','notGraded'};
        graded_unit_sum = cell(size(unit_h.PD));
        graded_unit_sum(unit_h.PD & unit_h.OD) = graded_name(1);
        graded_unit_sum(unit_h.PD & ~unit_h.OD) = graded_name(2);
        graded_unit_sum(~unit_h.PD & unit_h.OD) = graded_name(3);
        graded_unit_sum(~unit_h.PD & ~unit_h.OD) = graded_name(4);
        graded_unit_num = sum([unit_h.PD & unit_h.OD,unit_h.PD & ~unit_h.OD,~unit_h.PD & unit_h.OD,~unit_h.PD & ~unit_h.OD]);
    
        save([result_path new_filename{nf} '_graded_unit.mat'],'unit_estimate','graded_unit_num','graded_unit_sum','unit_h')
    else
        load([result_path new_filename{nf} '_graded_unit.mat'])
    end
    
    % add graded_code to spk_table
    tmp_graded_spk = spk_table.tune_code;
    tmp_graded = unitPrfDir.tune_code;
    tmp_maxDir = unitPrfDir.maxDir;
    unit_graded_code = unitPrfDir.tune_code;
    for uu = 1:nUnit_sel
        cur_unit = spk_table.unitName == unit_estimate.unitName(uu);
        sel_unit = unitPrfDir.unitName == unit_estimate.unitName(uu);
        if unit_h.PD(uu) && unit_h.OD(uu), tmp_code = 3;
        elseif unit_h.PD(uu) && ~unit_h.OD(uu), tmp_code = 2;
        elseif ~unit_h.PD(uu) && unit_h.OD(uu), tmp_code = 1;
        else tmp_code = 0;
        end
        tmp_graded_spk(cur_unit) = tmp_code;
        tmp_graded(sel_unit) = tmp_code;
        tmp_maxDir(uu) = unitPrfDir.maxDir(sel_unit);
        unit_graded_code(uu) = tmp_code;
    end
    spk_table.graded_code = tmp_graded_spk;
    unitPrfDir.graded_code = tmp_graded;
    unit_estimate.maxDir = tmp_maxDir;
    unit_estimate.graded_code = unit_graded_code;
    

    spk_mean = unique(outerjoin(spk_mean,spk_table,'Keys',{'unitName',biasName,setting.dtun},'Type','left','RightVariables',{'graded_code'}),'rows');

    fr_name = 'mean_firingRate';
    yLim = 30;
    tmp_suffix = '';
    
    unit_act = [];
    for uu = 1:nUnit_sel
        cur_unit = spk_mean.unitName == unitPrfDir{uu,1};
        cur_PD = round(unitPrfDir{uu,2});
        cur_maxDir = unique(spk_mean.maxDir(cur_unit));
        tmp_unit = spk_mean(cur_unit,{biasName,setting.dtun,fr_name});
        tmp_graded = spk_mean(cur_unit,{'graded_code'});
        tmpComponents = table2array(unstack(tmp_unit,fr_name,setting.dtun));
        maxInd = find(dir_list == cur_maxDir);
        tmp_circ = circshift((1:nDir)',1-maxInd);
        for dd = 1:nDir
            tmp_ind = tmp_circ(dd);
            tmp_act = array2table(tmpComponents(:,tmp_ind+1),'VariableNames',{'firingRate'});
            tmp_cond = array2table(cond_list,'VariableNames',{biasName});
            tmp_name = array2table(repmat(unitPrfDir{uu,1},size(tmp_act)),'VariableNames',{'unitName'});
            tmp_dir = array2table(repmat(relat_dir_name(dd),size(tmp_act)),'VariableNames',{'relatDir'});
            tmp_graded = tmp_graded(1:size(tmp_act),:);
            unit_act = [unit_act;[tmp_name,tmp_act,tmp_dir,tmp_cond,tmp_graded]];
        end
    end
      
    for gg = unique(spk_mean.graded_code)'
        tmp_act = unit_act(unit_act.graded_code == gg,:);
        unstack_unit_act = unstack(tmp_act,'firingRate','relatDir');
        spk_pop_mean = zeros(nDir,ncond);
        spk_pop_se = zeros(nDir,ncond);
        spk_h = zeros(nDir,ncond-1);
        spk_p = nan(nDir,ncond-1);
        test_norm1 = zeros(nDir,ncond);
        test_norm2 = nan(nDir,ncond);

        ii = 1;
        for dd = 1:length(relat_dir_name)
            tmp_unstack = unstack(unstack_unit_act(:,{'unitName',relat_dir_name{dd},'BiasTypeAbs'}),relat_dir_name{dd},'BiasTypeAbs');
            [spk_pop_mean(dd,:), ~,spk_pop_se(dd,:)] = imean(tmp_unstack{:,2:end},1);
            spk_pop_median(dd,:) = imedian(tmp_unstack{:,2:end},1);
            spk_pop_seM(dd,:) = spk_pop_se(dd,:)*sqrt(pi/2);
            
            for cc = 1:ncond-1
                [spk_p(dd,cc),spk_h(dd,cc)] = signrank(tmp_unstack{:,cc+1},tmp_unstack{:,cc+2});
            end
        end
        
        maxFir = max(max(spk_pop_median));
        maxFir = maxFir + 0.2*maxFir;
        h = figure('Name',[new_filename{nf} '_popAvg_errorBar' tmp_suffix],'Position',get(0,'ScreenSize'));
%                 errorbar(spk_pop_mean',spk_pop_se','LineWidth',2)
        errorbar(spk_pop_median',spk_pop_seM','LineWidth',2)
%                 for dd = 1:4
%                     errorbar(bias_strength,spk_pop_median(dd,:),spk_pop_median(dd,:)-spk_pop_lo(dd,:),spk_pop_up(dd,:)-spk_pop_median(dd,:),'LineWidth',2)
%                     hold on
%                 end
        text(1,maxFir*4/5,num2str(spk_p));
        ylim([0 maxFir]);
        legend(relat_dir_name)
        title({['FR of ' num2str(nUnit_sel) ' units at different bias levels ( ' setting.dtun ' rel. to prefDir/ ' tmp_suffix(2:end) ' )'];new_filename{nf}})
        
    saveas(h,[temp_fig_path new_filename{nf} '_popAvg_errorBar' tmp_suffix '.fig']);
    save2pdf([temp_fig_path new_filename{nf} '_popAvg_errorBar' tmp_suffix '.pdf'],h);
    caf
        
        maxFir = max(max(spk_pop_mean));
        maxFir = maxFir + 0.2*maxFir;
        expand_spk_pop_mean = [fliplr(spk_pop_mean(3:4,2:4)),spk_pop_mean(1:2,:)];
        expand_spk_pop_se = [fliplr(spk_pop_se(3:4,2:4)),spk_pop_se(1:2,:)];
        expand_spk_p = [fliplr(spk_p(3:4,:)),spk_p(1:2,:)];
        h = figure('Name',[new_filename{nf} '_popAvg_errorBar2' tmp_suffix],'Position',get(0,'ScreenSize'));
        errorbar(expand_spk_pop_mean',expand_spk_pop_se','LineWidth',2)
        text(0.05,maxFir*1.5/5,['p OD-OP = ' num2str(expand_spk_p(1,:))]);
        text(0.05,maxFir/5,['p Orth = ' num2str(expand_spk_p(2,:))]);
        text(-0.7,maxFir*4.5/5,['avg OD-PD = ' num2str(expand_spk_pop_mean(1,:))]);
        text(-0.7,maxFir*4/5,['avgO rth = ' num2str(expand_spk_pop_mean(2,:))]);
        ylim([0 maxFir]);
        legend({'OD-PD','Orth'})
        title({['FR of ' num2str(nUnit_sel) ' units at different bias levels ( ' setting.dtun ' rel. to prefDir/ ' tmp_suffix(2:end) ' )'];new_filename{nf}})
        
    saveas(h,[temp_fig_path new_filename{nf} '_popAvg_errorBar2' tmp_suffix '.fig']);
    save2pdf([temp_fig_path new_filename{nf} '_popAvg_errorBar2' tmp_suffix '.pdf'],h);
    caf

    end
    
    
    unit_estimate.session = cellfun(@(x) x(1:10),cellstr(unit_estimate.unitName),'UniformOutput',0);
    
    g = gramm('x',unit_estimate.OD,'y',unit_estimate.PD,'color',unit_estimate.session);
    g.set_names('color','Session')
    g.geom_point()
    g.stat_glm('geom','line','disp_fit','true')
%     g.stat_ellipse('type','95percentile','geom','area','patch_opts',{'FaceAlpha',0.1,'LineWidth',0.5})
    g.axe_property('DataAspectRatio',[1 1 1])
    g.geom_hline('yintercept',0,'style','k:')
    g.geom_vline('xintercept',0,'style','k:')
    g.set_names('x','OD slope','y', 'PD slope')
    g.set_title({['graded activities of ' num2str(graded_unit_num) ' units'];new_filename{nf}})
    h = figure('Name',[new_filename{nf} '_singleUnits_gradedTypeBySession'],'Position',get(0,'ScreenSize'));
    g.draw()
    saveas(h,[temp_fig_path new_filename{nf} '_singleUnits_gradedTypeBySession.fig']);
    save2pdf([temp_fig_path new_filename{nf} '_singleUnits_gradedTypeBySession.pdf'],h);
    caf
    
    
    g = gramm('x',unit_estimate.OD,'y',unit_estimate.PD,'color',unit_estimate.maxDir);
    g.set_names('color','maxDir')
    g.geom_point()
%     g.stat_glm('geom','line','disp_fit','true')
    g.stat_ellipse('type','95percentile','geom','area','patch_opts',{'FaceAlpha',0.1,'LineWidth',0.5})
    g.axe_property('DataAspectRatio',[1 1 1])
    g.geom_hline('yintercept',0,'style','k:')
    g.geom_vline('xintercept',0,'style','k:')
    g.set_names('x','OD slope','y', 'PD slope')
    g.set_title({['graded activities of ' num2str(graded_unit_num) ' units'];new_filename{nf}})
    h = figure('Name',[new_filename{nf} '_singleUnits_gradedTypeByMaxDir'],'Position',get(0,'ScreenSize'));
    g.draw()
    saveas(h,[temp_fig_path new_filename{nf} '_singleUnits_gradedTypeByMaxDir.fig']);
    save2pdf([temp_fig_path new_filename{nf} '_singleUnits_gradedTypeByMaxDir.pdf'],h);
    caf
    
    unit_estimate.electrode = cellfun(@(x) x(1:end-1), cellstr(unit_estimate.unitName),'UniformOutput',0);
    g = gramm('x',unit_estimate.OD,'y',unit_estimate.PD,'color',unit_estimate.electrode);
    g.geom_point()
%     g.stat_glm('geom','line','disp_fit','true')
    %     g.stat_ellipse('type','95percentile','geom','area','patch_opts',{'FaceAlpha',0.1,'LineWidth',0.5})
    g.axe_property('DataAspectRatio',[1 1 1])
    g.geom_hline('yintercept',0,'style','k:')
    g.geom_vline('xintercept',0,'style','k:')
    g.set_names('x','OD slope','y', 'PD slope')
    g.set_title({['graded activities of ' num2str(graded_unit_num) ' units'];new_filename{nf}})
    h = figure('Name',[new_filename{nf} '_singleUnits_gradedTypeByElectrode'],'Position',get(0,'ScreenSize'));
    g.draw()
    saveas(h,[temp_fig_path new_filename{nf} '_singleUnits_gradedTypeByElectrode.fig']);
    save2pdf([temp_fig_path new_filename{nf} '_singleUnits_gradedTypeByElectrode.pdf'],h);
    
    g = gramm('x',unit_estimate.OD,'y',unit_estimate.PD,'color',graded_unit_sum);
    g.set_names('color','gradedType')
    g.geom_point()
%     g.stat_glm('geom','line','disp_fit','true')
    g.stat_ellipse('type','95percentile','geom','area','patch_opts',{'FaceAlpha',0.1,'LineWidth',0.5})
    g.axe_property('DataAspectRatio',[1 1 1])
    g.geom_hline('yintercept',0,'style','k:')
    g.geom_vline('xintercept',0,'style','k:')
    g.set_names('x','OD slope','y', 'PD slope')
    g.set_title({['graded activities of ' num2str(graded_unit_num) ' units'];new_filename{nf}})
    h = figure('Name',[new_filename{nf} '_singleUnits_gradedType'],'Position',get(0,'ScreenSize'));
    g.draw()
    saveas(h,[temp_fig_path new_filename{nf} '_singleUnits_gradedType.fig']);
    save2pdf([temp_fig_path new_filename{nf} '_singleUnits_gradedType.pdf'],h);
    caf
     
    [unit_estimate.theta,unit_estimate.rho] = cart2pol(unit_estimate.OD,unit_estimate.PD);
    unit_estimate.deg = rad2deg(unit_estimate.theta);
    unit_estimate.deg(unit_estimate.deg<0) = 360 + unit_estimate.deg(unit_estimate.deg<0);
    nboot = 10000;
    [dip, p] = HartigansDipSignifTest(unit_estimate.deg(unit_estimate.deg > 45 & unit_estimate.deg < 255), nboot);
   
    h = figure('Name',[new_filename{nf} '_singleUnits_graded_hist'],'Position',get(0,'ScreenSize'));
    hist(unit_estimate.deg(unit_estimate.deg > 45 & unit_estimate.deg < 255),45:5:255)
    title([new_filename{nf} '_singleUnits_graded | dip=',num2str(dip,3), ', p=',num2str(p,3)],'Interpreter','none')
    set(gca,'XTick',90:90:180,'XTickLabel',{'PD+','OD-'})
    ylim([0 40])
    line([90 90],[0 40],'LineStyle',':')
    line([180 180],[0 40],'LineStyle',':')
    saveas(h,[temp_fig_path new_filename{nf} '_singleUnits_graded_hist.fig']);
    save2pdf([temp_fig_path new_filename{nf} '_singleUnits_graded_hist.pdf'],h);
    caf
    end