classdef PKpredictor
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name                = '';
        gcn                 = '';
        smiles              = '';
        pKa                 = struct('acidic',{},'basic',{},'cmpd_form',{},'source',{},'study',{},'lnb_ref',{},'date_experiment',{},'date_record',{},'date_updated',{},'comments',{}); % pKas
        logP                = struct('cmpd',{},'value',{},'source',{}); % logP;
        bp                  = struct('cmpd',{},'value',{},'species',{},'details',{},'source',{},'flag',{},'study',{},'lnb_ref',{},'date_experiment',{},'date_record',{},'date_updated',{},'comments',{}) % Blood:Plasma Ratio;
        fup                 = struct('cmpd',{},'value',{},'species',{},'details',{},'source',{},'method',{},'study',{},'lnb_ref',{},'date_experiment',{},'date_record',{},'date_updated',{},'comments',{});          % Fraction unbound in plasma;
        bind                = struct('cmpd',{},'matrix',{},'fu',{},'species',{},'details',{},'source',{},'method',{},'study',{},'lnb_ref',{},'date_experiment',{},'date_record',{},'date_updated',{},'comments',{}); % Protein Binding;
        met                 = struct('cmpd',{},'Modifier',{},'CL',{},'unit',{},'species',{},'system',{},'details',{},'source',{},'study',{},'lnb_ref',{},'date_experiment',{},'date_record',{},'date_updated',{},'comments',{});% In vitro metabolism;
        pk_kate             = struct('cmpd',{},'species',{},'strain',{},'sex',{},'n',{},'dose_mgkg',{},'route',{},'formulation',{},'matrix',{},'F',{},'F_std',{},'Cmax_ngml',{},'Cmax_sd',{},'Tmax_hr',{},'Tmax_sd',{},'AUC_Inf_nghrmL',{},'AUC_Inf_std',{},'TBdy_CL_mlmin_kg',{},'TBdy_CL_std',{},'Vdss_Lkg',{},'Vdss_std',{},'T1_2_hr',{},'T1_2_std',{},'project',{},'study',{},'study_design',{},'pk_start_date',{},'pk_date_record',{},'pk_date_updated',{},'comments',{});% In vivo PK data from KATE;
        cmpd_inputs         = struct('cmpd',{},'logP',{},'logP_source',{},'logD',{},'logD_pH',{},'logD_source',{},'Acidic_pKa',{},'Acidic_pKa_source',{},'Basic_pKa',{},'Basic_pKa_source',{});
        species_inputs      = struct('cmpd',{},'species',{},'CLint_mic',{},'CLint_mic_pctCV',{},'CLint_mic_source',{},'CLint_hep',{},'CLint_hep_pctCV',{},'CLint_hep_source',{},'fup',{},'fup_pctCV',{},'fup_source',{},'bp',{},'bp_pctCV',{},'bp_source',{},'fub',{},'fub_pctCV',{},'fub_source',{},'fuinc_hep',{},'fuinc_hep_source',{},'fuinc_mic',{},'fuinc_mic_source',{});
        ivive               = struct('cmpd',{},'species',{},'type',{},'system',{},'CLp_mLminkg',{},'CLb_mLminkg',{});
        ivive_results       = struct('cmpd',{},'species',{},'strain',{},'type',{},'system',{},'matrix',{},'CL_mLminkg',{},'CLref_mLminkg',{},'fold_error',{},'CLref_pctCV',{},'CLref_n',{},'dose_mgkg',{},'study',{},'study_design',{});
        lbf                 = struct('cmpd',{},'species',{},'strain',{},'sex',{},'matrix',{},'CL_mLminkg',{},'Thalf_hr',{},'type',{},'HuCLp_mLminkg',{},'HuCLb_mLminkg',{},'HuThalf_hr',{},'animal_id',{},'study',{},'study_design',{});
        single_allometry    = struct('cmpd',{},'species',{},'strain',{},'sex',{},'matrix',{},'CLb_mLminkg',{},'CLp_mLminkg',{},'Vss_b_Lkg',{},'Vss_p_Lkg',{},'type',{},'HuCLp_mLminkg',{},'HuCLb_mLminkg',{},'HuVss_b_Lkg',{},'HuVss_p_Lkg',{},'HuThalf_hr',{},'animal_id',{},'study',{},'study_design',{});
        allometry_data      = struct('cmpd',{},'species',{},'BW_kg',{},'BrW_g',{},'MLP_yrs',{},'matrix',{},'CL_mLminkg',{},'Vss_Lkg',{},'CLb_mLmin',{},'CLp_mLmin',{},'CLbU_mLmin',{},'CLpU_mLmin',{},'Vss_b_L',{},'Vss_p_L',{},'VssU_b_L',{},'VssU_p_L',{},'CLb_x_MLP',{},'CLp_x_MLP',{},'CLb_x_BrW',{},'CLp_x_BrW',{},'CLbU_x_MLP',{},'CLpU_x_MLP',{},'CLbU_x_BrW',{},'CLpU_x_BrW',{});
        multi_allometry     = struct('cmpd',{},'parameter',{},'method',{},'type',{},'value_blood',{},'value_plasma',{},'R2_fit',{},'a',{},'b',{},'equation',{},'preferred_cl_res',{},'preferred_cl_nonres',{});
        phys                = struct('species',{},'BW_kg',{},'LBF_mLminkg',{},'SF_gLiver_kgBW',{},'GFR_mLminkg',{},'BrW_g',{},'LiverW_g',{},'MLP_yrs',{});
        calc_opt            = struct('parameter_operation',{},'pk_data_operation',{},'use_iv_data',{});
    end
    
    methods
        function obj = PKpredictor(gcn,varargin)
            if isa(gcn,'smtb.DrugProps')
                dp = gcn;
            elseif isa(gcn,'char')
                dp = smtb.DrugProps(gcn);
            else
                error('Input gcn must be either a gsk compound number as char or DrugProps object.')
            end
            
            obj.name = dp.name;
            obj.gcn = dp.gcn;
            obj.smiles = dp.smiles;
            
            obj.pKa = dp.pKa;
            obj.logP = dp.logP;
            obj.bp = dp.bp;
            obj.fup = dp.fup;
            obj.bind = dp.bind;
            obj.met = dp.met;
            obj.pk_kate = dp.pk_kate;
            
            species = {'Human','Minipig','Dog','Rat','Mouse','Monkey'};
            
            cmpd_opts = struct('logP',[],'logP_source','default',...
                'logD',[],'logD_pH',7.4,'logD_source','calculated',...
                'Acidic_pKa',[],'Acidic_pKa_source','default',...
                'Basic_pKa',[],'Basic_pKa_source','default');
            
            sp_opts = struct('species',species,...
                'fup',[],'adjusted_fup',[],'fup_type','experimental','fup_pctCV',[],'fup_source','default',...
                'bp',[],'bp_pctCV',[],'bp_source','default',...
                'fub',[],'fub_pctCV',[],'fub_source','default',...
                'CLint_mic',[],'CLint_mic_pctCV',[],'CLint_mic_source','default',...
                'CLint_hep',[],'CLint_hep_pctCV',[],'CLint_hep_source','default',...
                'CLint_S9',[],'CLint_S9_pctCV',[],'CLint_S9_source','default',...
                'fuinc_mic',1,'fuinc_mic_source','default',...
                'fuinc_hep',1,'fuinc_hep_source','default',...
                'fuinc_S9',1,'fuinc_S9_source','default');
            
            calc_opts = struct('iv_data','all',...
                'param_op','arithmetic mean',...
                'pk_data_op','geometric mean');
            
            phys = [];
            
            for i = 1:2:length(varargin)
                switch strrep(strrep(lower(varargin{i}),'_',''),' ','')
                    case {'phys','physiology'}
                        phys = smtb.useful.typecheck(varargin{i+1},'struct','Physiology');
                    case 'ivdata'
                        calc_opts.iv_data = smtb.useful.typecheck(varargin{i+1},'char','IV data');
                        switch strrep(strrep(lower(varargin{i+1}),'_',''),' ','')
                            case {'discrete','cassette','all'}
                                calc_opts.iv_data = smtb.useful.typecheck(varargin{i+1},'char','IV data');
                            otherwise
                                error('IV data choice not valid. Valid options are discrete, cassette, and median')
                        end
                    case {'params','parameters','parameteroperation'}
                        calc_opts.param_op = smtb.useful.typecheck(varargin{i+1},'char','parameters');
                        switch strrep(strrep(lower(varargin{i+1}),'_',''),' ','')
                            case 'arithmeticmean'
                                calc_opts.param_op = 'arithmetic mean';
                            case 'geometricmean'
                                calc_opts.param_op = 'geometric mean';
                            case 'median'
                                calc_opts.param_op = 'median';
                            otherwise
                                error('Parameter operation input not valid. Valid options are arithmetic mean, geometric mean, and median')
                        end
                    case {'pkdata','pkdataoperation'}
                        calc_opts.pk_data_op = smtb.useful.typecheck(varargin{i+1},'char','PK data operation');
                        switch strrep(strrep(lower(varargin{i+1}),'_',''),' ','')
                            case 'arithmeticmean'
                                calc_opts.pk_data_op = 'arithmetic mean';
                            case 'geometricmean'
                                calc_opts.pk_data_op = 'geometric mean';
                            case 'median'
                                calc_opts.pk_data_op = 'median';
                            otherwise
                                error('PK data operation input not valid. Valid options are arithmetic mean, geometric mean, and median')
                        end
                    case 'logp'
                        if isnumeric(varargin{i+1})
                            cmpd_opts.logP = smtb.useful.typecheck(varargin{i+1},'double','logP');
                            cmpd_opts.logP_source = 'User Defined';
                        elseif contains(varargin{i+1},'logd','IgnoreCase',true)
                            cmpd_opts.logD_source = smtb.useful.typecheck(varargin{i+1},'char','logD source');
                            cmpd_opts.logP_source = 'calculated';
                        else
                            cmpd_opts.logP_source = smtb.useful.typecheck(varargin{i+1},'char','logP source');
                        end
                    case 'logd'
                        if isnumeric(varargin{i+1})
                            cmpd_opts.logD = smtb.useful.typecheck(varargin{i+1},'double','logD');
                            cmpd_opts.logD_source = 'User Defined';
                            cmpd_opts.logP_source = 'calculated';
                        else
                            cmpd_opts.logD_source = smtb.useful.typecheck(varargin{i+1},'char','logD source');
                            cmpd_opts.logP_source = 'calculated';
                        end
                    case 'logdph'
                        cmpd_opts.logD_pH = smtb.useful.typecheck(varargin{i+1},'double','logD pH');
                    case 'acidicpka'
                        if isnumeric(varargin{i+1})
                            cmpd_opts.Acidic_pKa = smtb.useful.typecheck(varargin{i+1},'double','Acidic pKa');
                            cmpd_opts.Acidic_pKa_source = 'User Defined';
                        else
                            cmpd_opts.Acidic_pKa_source = smtb.useful.typecheck(varargin{i+1},'char','Acidic pKa source');
                        end
                    case 'basicpka'
                        if isnumeric(varargin{i+1})
                            cmpd_opts.Basic_pKa = smtb.useful.typecheck(varargin{i+1},'double','Basic pKa');
                            cmpd_opts.Basic_pKa_source = 'User Defined';
                        else
                            cmpd_opts.Basic_pKa_source = smtb.useful.typecheck(varargin{i+1},'char','Basic pKa source');
                        end
                    otherwise
                        % Parse species inputs
                        if ~isempty(sscanf(lower(varargin{i}),'fup_%s'))
                            sp      = sscanf(lower(varargin{i}),'fup_%s');
                            [~,ia]  = intersect(lower({sp_opts.species}),lower(sp));
                            if isempty(ia), error('Unrecognized species "%s" in "%s"',sp,varargin{i})
                            else
                                if isnumeric(varargin{i+1})
                                    sp_opts(ia).fup    = smtb.useful.typecheck(varargin{i+1},'double',varargin{i+1});
                                    sp_opts(ia).fup_source = 'User Defined';
                                else
                                    sp_opts(ia).fup_source = smtb.useful.typecheck(varargin{i+1},'char',varargin{i+1});
                                end
                            end
                        elseif ~isempty(sscanf(lower(varargin{i}),'bp_%s'))
                            sp = sscanf(lower(varargin{i}),'bp_%s');
                            [~,ia]  = intersect(lower({sp_opts.species}),lower(sp));
                            if isempty(ia), error('Unrecognized species "%s" in "%s"',sp,varargin{i})
                            else
                                if isnumeric(varargin{i+1})
                                    sp_opts(ia).bp    = smtb.useful.typecheck(varargin{i+1},'double',varargin{i+1});
                                    sp_opts(ia).bp_source = 'User Defined';
                                else
                                    sp_opts(ia).bp_source = smtb.useful.typecheck(varargin{i+1},'char',varargin{i+1});
                                end
                            end
                        elseif ~isempty(sscanf(lower(varargin{i}),'fub_%s'))
                            sp = sscanf(lower(varargin{i}),'fub_%s');
                            [~,ia]  = intersect(lower({sp_opts.species}),lower(sp));
                            if isempty(ia), error('Unrecognized species "%s" in "%s"',sp,varargin{i})
                            else
                                if isnumeric(varargin{i+1})
                                    sp_opts(ia).fub    = smtb.useful.typecheck(varargin{i+1},'double',varargin{i+1});
                                    sp_opts(ia).fub_source = 'User Defined';
                                else
                                    sp_opts(ia).fub_source = smtb.useful.typecheck(varargin{i+1},'char',varargin{i+1});
                                end
                            end
                        elseif ~isempty(sscanf(lower(varargin{i}),'clint_mic_%s'))
                            sp = sscanf(lower(varargin{i}),'clint_mic_%s');
                            [~,ia]  = intersect(lower({sp_opts.species}),lower(sp));
                            if isempty(ia), error('Unrecognized species "%s" in "%s"',sp,varargin{i})
                            else
                                if isnumeric(varargin{i+1})
                                    sp_opts(ia).CLint_mic    = smtb.useful.typecheck(varargin{i+1},'double',varargin{i+1});
                                    sp_opts(ia).CLint_mic_source = 'User Defined';
                                else
                                    sp_opts(ia).CLint_mic_source = smtb.useful.typecheck(varargin{i+1},'char',varargin{i+1});
                                end
                            end
                        elseif ~isempty(sscanf(lower(varargin{i}),'fuinc_mic_%s'))
                            sp = sscanf(lower(varargin{i}),'fuinc_mic_%s');
                            [~,ia]  = intersect(lower({sp_opts.species}),lower(sp));
                            if isempty(ia), error('Unrecognized species "%s" in "%s"',sp,varargin{i})
                            else
                                if isnumeric(varargin{i+1})
                                    sp_opts(ia).fuinc_mic    = smtb.useful.typecheck(varargin{i+1},'double',varargin{i+1});
                                    sp_opts(ia).fuinc_mic_source = 'User Defined';
                                else
                                    sp_opts(ia).fuinc_mic_source = smtb.useful.typecheck(varargin{i+1},'char',varargin{i+1});
                                end
                            end
                        elseif ~isempty(sscanf(lower(varargin{i}),'clint_hep_%s'))
                            sp = sscanf(lower(varargin{i}),'clint_hep_%s');
                            [~,ia]  = intersect(lower({sp_opts.species}),lower(sp));
                            if isempty(ia), error('Unrecognized species "%s" in "%s"',sp,varargin{i})
                            else
                                if isnumeric(varargin{i+1})
                                    sp_opts(ia).CLint_hep    = smtb.useful.typecheck(varargin{i+1},'double',varargin{i+1});
                                    sp_opts(ia).CLint_hep_source = 'User Defined';
                                else
                                    sp_opts(ia).CLint_hep_source = smtb.useful.typecheck(varargin{i+1},'char',varargin{i+1});
                                end
                            end
                        elseif ~isempty(sscanf(lower(varargin{i}),'fuinc_hep_%s'))
                            sp = sscanf(lower(varargin{i}),'fuinc_hep_%s');
                            [~,ia]  = intersect(lower({sp_opts.species}),lower(sp));
                            if isempty(ia), error('Unrecognized species "%s" in "%s"',sp,varargin{i})
                            else
                                if isnumeric(varargin{i+1})
                                    sp_opts(ia).fuinc_hep    = smtb.useful.typecheck(varargin{i+1},'double',varargin{i+1});
                                    sp_opts(ia).fuinc_hep_source = 'User Defined';
                                else
                                    sp_opts(ia).fuinc_hep_source = smtb.useful.typecheck(varargin{i+1},'char',varargin{i+1});
                                end
                            end
                        elseif ~isempty(sscanf(lower(varargin{i}),'fuinc_S9_%s'))
                            sp = sscanf(lower(varargin{i}),'fuinc_S9_%s');
                            [~,ia]  = intersect(lower({sp_opts.species}),lower(sp));
                            if isempty(ia), error('Unrecognized species "%s" in "%s"',sp,varargin{i})
                            else
                                if isnumeric(varargin{i+1})
                                    sp_opts(ia).fuinc_S9    = smtb.useful.typecheck(varargin{i+1},'double',varargin{i+1});
                                    sp_opts(ia).fuinc_S9_source = 'User Defined';
                                else
                                    sp_opts(ia).fuinc_S9_source = smtb.useful.typecheck(varargin{i+1},'char',varargin{i+1});
                                end
                            end
                            % following lines allow one fuinc to be defined
                            % defined for both hep and mic
                        elseif ~isempty(sscanf(lower(varargin{i}),'fuinc_%s'))
                            sp = sscanf(lower(varargin{i}),'fuinc_%s');
                            [~,ia]  = intersect(lower({sp_opts.species}),lower(sp));
                            if isempty(ia), error('Unrecognized species "%s" in "%s"',sp,varargin{i})
                            else
                                if isnumeric(varargin{i+1})
                                    sp_opts(ia).fuinc_hep    = smtb.useful.typecheck(varargin{i+1},'double',varargin{i+1});
                                    sp_opts(ia).fuinc_mic    = smtb.useful.typecheck(varargin{i+1},'double',varargin{i+1});
                                    sp_opts(ia).fuinc_hep_source = 'User Defined';
                                    sp_opts(ia).fuinc_mic_source = 'User Defined';
                                else
                                    sp_opts(ia).fuinc_hep_source = smtb.useful.typecheck(varargin{i+1},'char',varargin{i+1});
                                    sp_opts(ia).fuinc_mic_source = smtb.useful.typecheck(varargin{i+1},'char',varargin{i+1});
                                end
                            end
                        end
                end
            end
            
            %%%% Define compound inputs %%%%
            % logP
            if ~(strcmpi(cmpd_opts.logP_source,'User Defined') || strcmpi(cmpd_opts.logP_source,'calculated'))
                [logP,logP_source] = getPreferredValue(dp,'logP',cmpd_opts.logP_source);
                if contains(logP_source,'logd','IgnoreCase',true)
                    cmpd_opts.logD_source = logP_source;
                    cmpd_opts.logP_source = 'calculated';
                else
                    cmpd_opts.logP = logP;
                    cmpd_opts.logP_source = logP_source;
                end
            end
            % logD
            if ~(strcmpi(cmpd_opts.logD_source,'User Defined') || strcmpi(cmpd_opts.logD_source,'calculated'))
                [logD,logD_source] = getPreferredValue(dp,'logP',cmpd_opts.logD_source);
                if contains(logD_source,'logP','IgnoreCase',true) || contains(logD_source,'QSAR Biobyte/Daylight','IgnoreCase',true)
                    cmpd_opts.logP = logD;
                    cmpd_opts.logP_source = logD_source;
                    cmpd_opts.logD_source = 'calculated';
                else
                    cmpd_opts.logD = logD;
                    cmpd_opts.logD_source = logD_source;
                end
            end
            
            % Acidic pKa
            if ~strcmpi(cmpd_opts.Acidic_pKa_source,'User Defined')
                [pKa_a,pKa_a_source] = getPreferredValue(dp,'Acidic_pKa',cmpd_opts.Acidic_pKa_source);
                cmpd_opts.Acidic_pKa = pKa_a;
                cmpd_opts.Acidic_pKa_source = pKa_a_source;
            end
            
            % Basic pKa
            if ~strcmpi(cmpd_opts.Basic_pKa_source,'User Defined')
                [pKa_b,pKa_b_source] = getPreferredValue(dp,'Basic_pKa',cmpd_opts.Basic_pKa_source);
                cmpd_opts.Basic_pKa = pKa_b;
                cmpd_opts.Basic_pKa_source = pKa_b_source;
            end
            
            % calculate logD/P
            if strcmpi(cmpd_opts.logD_source,'calculated')
                cmpd_opts.logD = smtb.pbpk.logDcalc(cmpd_opts.logP,cmpd_opts.logD_pH,cmpd_opts.Acidic_pKa,cmpd_opts.Basic_pKa);
            elseif strcmpi(cmpd_opts.logP_source,'calculated')
                cmpd_opts.logP = smtb.pbpk.logPcalc(cmpd_opts.logD,cmpd_opts.logD_pH,cmpd_opts.Acidic_pKa,cmpd_opts.Basic_pKa);
            else
                error('logP or logD need to be source "calculated"')
            end
            
            % add all compound properties to object
            obj = obj.addPropertyValue('cmpd_inputs','logP',cmpd_opts.logP,'logP_source',cmpd_opts.logP_source,...
                'logD',cmpd_opts.logD,'logD_pH',cmpd_opts.logD_pH,'logD_source',cmpd_opts.logD_source,...
                'Acidic_pKa',cmpd_opts.Acidic_pKa,'Acidic_pKa_source',cmpd_opts.Acidic_pKa_source,...
                'Basic_pKa',cmpd_opts.Basic_pKa,'Basic_pKa_source',cmpd_opts.Basic_pKa_source);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%% Define species-specific inputs %%%%
            for i = 1:length(species)
                idx = strcmpi({sp_opts.species},species{i});
                % fup
                if ~strcmpi(sp_opts(idx).fup_source,'User Defined')
                    [value,fup_source,sp,~,~,~,~,fup_vals] = getPreferredValue(dp,'fup',species{i},sp_opts(idx).fup_source,calc_opts.param_op);
                    sp_opts(idx).fup = value;
                    sp_opts(idx).fup_pctCV = (std(fup_vals)./mean(fup_vals)) * 100;
                    if contains(fup_source,'QSAR','IgnoreCase',true)
                        sp_opts(idx).fup_source = fup_source;
                    else
                        sp_opts(idx).fup_source = [fup_source ' ' upper(sp(1)) lower(sp(2:end))];
                    end
                end
                % bp
                switch lower(sp_opts(idx).bp_source)
                    case 'user defined' % do nothing
                    case 'fup/fub' % first check for fub values. if found, define fub here
                        if ~isempty(sp_opts(idx).fub) % if it is not empty, this means it was defined by user
                            sp_opts(idx).bp = sp_opts(idx).fup/sp_opts(idx).fub;
                        else
                            [value,fub_source,sp,~,~,~,~,fub_vals] = getPreferredValue(dp,'fub',species{i},sp_opts(idx).fub_source,calc_opts.param_op);
                            if isempty(fub_vals)
                                [value,bp_source,sp,~,~,~,~,bp_vals] = getPreferredValue(dp,'bp',species{i},sp_opts(idx).bp_source,calc_opts.param_op);
                                sp_opts(idx).bp = value;
                                sp_opts(idx).bp_pctCV = (std(bp_vals)./mean(bp_vals)) * 100;
                                if contains(bp_source,'QSAR','IgnoreCase',true)
                                    sp_opts(idx).bp_source = bp_source;
                                else
                                    sp_opts(idx).bp_source = [bp_source ' ' upper(sp(1)) lower(sp(2:end))];
                                end
                            else
                                sp_opts(idx).fub = value;
                                sp_opts(idx).bp = sp_opts(idx).fup./sp_opts(idx).fub;
                                sp_opts(idx).bp_source = 'fup/fub';
                                if contains(fub_source,'QSAR','IgnoreCase',true)
                                    sp_opts(idx).fub_source = fub_source;
                                else
                                    sp_opts(idx).fub_source = [fub_source ' ' upper(sp(1)) lower(sp(2:end))];
                                end
                            end
                        end
                    otherwise
                        [value,bp_source,sp,~,~,~,~,bp_vals] = getPreferredValue(dp,'bp',species{i},sp_opts(idx).bp_source,calc_opts.param_op);
                        sp_opts(idx).bp = value;
                        sp_opts(idx).bp_pctCV = (std(bp_vals)./mean(bp_vals)) * 100;
                        if contains(bp_source,'QSAR','IgnoreCase',true)
                            sp_opts(idx).bp_source = bp_source;
                        else
                            sp_opts(idx).bp_source = [bp_source ' ' upper(sp(1)) lower(sp(2:end))];
                        end
                end
                % fub
                if isempty(sp_opts(idx).fub)
                    switch lower(sp_opts(idx).fub_source)
                        case 'user defined' % do nothing
                        case 'fup/bp'
                            sp_opts(idx).fub = sp_opts(idx).fup./sp_opts(idx).bp;
                            sp_opts(idx).fub_source = 'fup/BP';
                        otherwise
                            [value,fub_source,sp,~,~,~,~,fub_vals] = getPreferredValue(dp,'fub',species{i},sp_opts(idx).fub_source,calc_opts.param_op);
                            sp_opts(idx).fub = value;
                            if isempty(fub_vals)
                                sp_opts(idx).fub = sp_opts(idx).fup./sp_opts(idx).bp;
                                sp_opts(idx).fub_source = 'fup/BP';
                            else
                                sp_opts(idx).fub_pctCV = (std(fub_vals)./mean(fub_vals)) * 100;
                                if contains(fub_source,'QSAR','IgnoreCase',true)
                                    sp_opts(idx).fub_source = fub_source;
                                else
                                    sp_opts(idx).fub_source = [fub_source ' ' upper(sp(1)) lower(sp(2:end))];
                                end
                            end
                    end
                end
                % CLint microsomes
                if ~strcmpi(sp_opts(idx).CLint_mic_source,'User Defined')
                    [value,cl_mic_source,sp,~,~,~,system,cl_mic_vals] = getPreferredValue(dp,'clint',species{i},{sp_opts(idx).CLint_mic_source,'KATE microsomes','QSAR ADMET','QSAR Met_Intrinsic_Clearance','QSAR Rat Intrinsic Clearance Class'},calc_opts.param_op);
                    if ~strcmpi(system,'microsomes')
                        [value,cl_mic_source,sp,~,~,~,system,cl_mic_vals] = getPreferredValue(dp,'clint','human',{sp_opts(idx).CLint_mic_source,'KATE microsomes','QSAR ADMET','QSAR Met_Intrinsic_Clearance','QSAR Rat Intrinsic Clearance Class'},calc_opts.param_op);
                    end
                    sp_opts(idx).CLint_mic = value;
                    if strcmpi(system,'microsomes')
                        sp_opts(idx).CLint_mic_pctCV = (std(cl_mic_vals)./mean(cl_mic_vals)) * 100;
                        if contains(cl_mic_source,'QSAR','IgnoreCase',true)
                            sp_opts(idx).CLint_mic_source = cl_mic_source;
                        else
                            sp_opts(idx).CLint_mic_source = [upper(system(1)) lower(system(2:end)) ' ' upper(sp(1)) lower(sp(2:end)) ' ' '(' cl_mic_source ')'];
                        end
                    end
                end
                % CLint hepatocytes
                if ~strcmpi(sp_opts(idx).CLint_hep_source,'User Defined')
                    [value,cl_hep_source,sp,~,~,~,system,cl_hep_vals] = getPreferredValue(dp,'clint',species{i},{sp_opts(idx).CLint_hep_source,'KATE hepatocytes'},calc_opts.param_op);
                    if ~strcmpi(system,'hepatocytes')
                        [value,cl_hep_source,sp,~,~,~,system,cl_hep_vals] = getPreferredValue(dp,'clint','human',{sp_opts(idx).CLint_hep_source,'KATE hepatocytes'},calc_opts.param_op);
                    end
                    if strcmpi(system,'hepatocytes')
                        sp_opts(idx).CLint_hep = value;
                        sp_opts(idx).CLint_hep_pctCV = (std(cl_hep_vals)./mean(cl_hep_vals)) * 100;
                        if contains(cl_hep_source,'QSAR','IgnoreCase',true)
                            sp_opts(idx).CLint_hep_source = cl_hep_source;
                        else
                            sp_opts(idx).CLint_hep_source = [upper(system(1)) lower(system(2:end)) ' ' upper(sp(1)) lower(sp(2:end)) ' ' '(' cl_hep_source ')'];
                        end
                    end
                end
                % CLint S9
                if ~strcmpi(sp_opts(idx).CLint_S9_source,'User Defined')
                    [value,cl_S9_source,sp,~,~,~,system,cl_S9_vals] = getPreferredValue(dp,'clint',species{i},{sp_opts(idx).CLint_hep_source,'KATE S9'},calc_opts.param_op);
                    if ~strcmpi(system,'s9')
                        [value,cl_S9_source,sp,~,~,~,system,cl_S9_vals] = getPreferredValue(dp,'clint','human',{sp_opts(idx).CLint_hep_source,'KATE S9'},calc_opts.param_op);
                    end
                    if strcmpi(system,'s9')
                        sp_opts(idx).CLint_S9 = value;
                        sp_opts(idx).CLint_S9_pctCV = (std(cl_S9_vals)./mean(cl_S9_vals)) * 100;
                        if contains(cl_hep_source,'QSAR','IgnoreCase',true)
                            sp_opts(idx).CLint_S9_source = cl_S9_source;
                        else
                            sp_opts(idx).CLint_S9_source = [upper(system(1)) lower(system(2:end)) ' ' upper(sp(1)) lower(sp(2:end)) ' ' '(' cl_S9_source ')'];
                        end
                    end
                end
                % fuinc microsomes
                if ~strcmpi(sp_opts(idx).fuinc_mic_source,'User Defined')
                    switch strrep(strrep(lower(sp_opts(idx).fuinc_mic_source),'_',''),' ','')
                        case 'kate'
                            [fuinc,source,sp,~,~,~,system] = getPreferredValue(dp,'fuinc_microsomes',species{i},'KATE',calc_opts.param_op);
                            if isempty(fuinc) || isnan(fuinc)
                                sp_opts(idx).fuinc_mic = 1;
                                sp_opts(idx).fuinc_mic_source = 'default';
                            else
                                sp_opts(idx).fuinc_mic = fuinc;
                                sp_opts(idx).fuinc_mic_source = [upper(system(1)) lower(system(2:end)) ' ' upper(sp(1)) lower(sp(2:end)) ' ' '(' source ')'];
                            end
                        case 'fup'
                            sp_opts(idx).fuinc_mic = sp_opts(idx).fup;
                            sp_opts(idx).fuinc_mic_source = 'fup';
                        case 'fub'
                            sp_opts(idx).fuinc_mic = sp_opts(idx).fub;
                            sp_opts(idx).fuinc_mic_source = 'fub';
                        case {'austinsimple','hallifax'}
                            sp_opts(idx).fuinc_mic = calcFuinc(cmpd_opts.logD,cmpd_opts.logD_pH,cmpd_opts.Acidic_pKa,cmpd_opts.Basic_pKa,sp_opts(idx).fuinc_mic_source);
                        case {'austin','austinmicrosomes'}
                            sp_opts(idx).fuinc_mic = calcFuinc(cmpd_opts.logD,cmpd_opts.logD_pH,cmpd_opts.Acidic_pKa,cmpd_opts.Basic_pKa,'Austin Microsomes');
                        case 'default' % do nothing
                        otherwise
                            error('fuinc microsome source not valid. Choose from Austin, Austin Simple, or Hallifax')
                    end
                end
                % fuinc hepatocytes
                if ~strcmpi(sp_opts(idx).fuinc_hep_source,'User Defined')
                    switch strrep(strrep(lower(sp_opts(idx).fuinc_hep_source),'_',''),' ','')
                        case 'kate'
                            [fuinc,source,sp,~,~,~,system] = getPreferredValue(dp,'fuinc_hepatocytes',species{i},'KATE',calc_opts.param_op);
                            if isempty(fuinc) || isnan(fuinc)
                                sp_opts(idx).fuinc_hep = 1;
                                sp_opts(idx).fuinc_hep_source = 'default';
                            else
                                sp_opts(idx).fuinc_hep = fuinc;
                                sp_opts(idx).fuinc_hep_source = [upper(system(1)) lower(system(2:end)) ' ' upper(sp(1)) lower(sp(2:end)) ' ' '(' source ')'];
                            end
                        case 'fup'
                            sp_opts(idx).fuinc_hep = sp_opts(idx).fup;
                            sp_opts(idx).fuinc_hep_source = 'fup';
                        case 'fub'
                            sp_opts(idx).fuinc_hep = sp_opts(idx).fub;
                            sp_opts(idx).fuinc_hep_source = 'fub';
                        case {'austin','austinhepatocytes'}
                            sp_opts(idx).fuinc_mic = calcFuinc(cmpd_opts.logD,cmpd_opts.logD_pH,cmpd_opts.Acidic_pKa,cmpd_opts.Basic_pKa,'Austin Hepatocytes');
                        case 'default' % do nothing
                        otherwise
                            error('fuinc microsome source not valid. Choose from Austin, Austin Simple, or Hallifax')
                    end
                end
                % add all values for species i
                obj = obj.addPropertyValue('species_inputs','species',species{i},...
                    'fup',sp_opts(idx).fup,'fup_pctCV',sp_opts(idx).fup_pctCV,'fup_source',sp_opts(idx).fup_source,...
                    'bp',sp_opts(idx).bp,'bp_pctCV',sp_opts(idx).bp_pctCV,'bp_source',sp_opts(idx).bp_source,...
                    'fub',sp_opts(idx).fub,'fub_pctCV',sp_opts(idx).fub_pctCV,'fub_source',sp_opts(idx).fub_source,...
                    'CLint_mic',sp_opts(idx).CLint_mic,'CLint_mic_pctCV',sp_opts(idx).CLint_mic_pctCV,'CLint_mic_source',sp_opts(idx).CLint_mic_source,...
                    'CLint_hep',sp_opts(idx).CLint_hep,'CLint_hep_pctCV',sp_opts(idx).CLint_hep_pctCV,'CLint_hep_source',sp_opts(idx).CLint_hep_source,...
                    'CLint_S9',sp_opts(idx).CLint_S9,'CLint_S9_pctCV',sp_opts(idx).CLint_S9_pctCV,'CLint_S9_source',sp_opts(idx).CLint_S9_source,...
                    'fuinc_mic',sp_opts(idx).fuinc_mic,'fuinc_mic_source',sp_opts(idx).fuinc_mic_source,...
                    'fuinc_hep',sp_opts(idx).fuinc_hep,'fuinc_hep_source',sp_opts(idx).fuinc_hep_source,...
                    'fuinc_S9',sp_opts(idx).fuinc_hep,'fuinc_S9_source',sp_opts(idx).fuinc_S9_source);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if isempty(phys)
                load('GastroplusPhysiology.mat','pbpk');
                phys = pbpk;
                clearvars pbpk
            end
            
            obj = obj.getPhysiologyConstants(obj,phys,species);
            
            obj = obj.calculateIVIVE(obj,species,calc_opts);
            
            obj = obj.LBFscaling(obj,calc_opts);
            
            obj = obj.SingleSpeciesAllometry(obj,calc_opts);
            
            obj = obj.BuildAllometryData(obj,calc_opts);
            
            obj = obj.MultiSpeciesAllometry(obj);
        end
        
        function obj = addPropertyValue(obj,prop,varargin)
            idx = length(obj.(prop))+1;
            for i = 1:2:length(varargin)
                obj.(prop)(idx).(varargin{i}) = varargin{i+1};
            end
        end
                
    end
    
    methods (Static,Hidden,Access=private)
        
        function obj = getPhysiologyConstants(obj,phys,species)
            for i = 1:length(species)
                physi = phys(strcmpi({phys.Species},species{i}));
                bw_kg = physi.Properties.Weight;
                tiss = physi.Tissues;
                % lbf units are stored as mL/s
                lbf = tiss(strcmpi({tiss.Name},'liver')).TissuePerfusion;
                % convert lbf units to mL/min/kg
                lbf = (lbf * 60) / bw_kg;
                % gfr units are stored as mL/s
                gfr = tiss(strcmpi({tiss.Name},'kidney')).GlomerularFiltrationRate;
                % convert gfr units to mL/min/kg
                gfr = (gfr * 60) / bw_kg;
                % liver weight (multiply density g/mL * volume mL)
                livw = tiss(strcmpi({tiss.Name},'liver')).Density * tiss(strcmpi({tiss.Name},'liver')).Volume;
                % scaling factor = g liver weight / kg bodyweight
                sf = livw / bw_kg;
                % brain weight (multiply density g/mL * volume mL)
                brw = tiss(strcmpi({tiss.Name},'brain')).Density * tiss(strcmpi({tiss.Name},'brain')).Volume;
                switch lower(species{i})
                    case 'human',       mlp = 93;
                    case 'mouse',       mlp = 2.7;
                    case 'rat',         mlp = 4.7;
                    case 'guinea pig',  mlp = 9;
                    case 'marmoset',    mlp = 10;
                    case 'dog',         mlp = 20;
                    case 'monkey',      mlp = 22;
                    case 'minipig',     mlp = 27;
                    case 'rabbit',      mlp = 8;
                end
                obj = obj.addPropertyValue('phys','species',species{i},'BW_kg',bw_kg,'LBF_mLminkg',lbf,'SF_gLiver_kgBW',sf,'GFR_mLminkg',gfr,'BrW_g',brw,'LiverW_g',livw,'MLP_yrs',mlp);
            end
        end
        
        function obj = calculateIVIVE(obj,species,calc_opts)
            %%%% calculate clearance from IVIVE %%%%
            for i = 1:length(species)
                sp_in = obj.species_inputs(strcmpi({obj.species_inputs.species},species{i}));
                sp_phys = obj.phys(strcmpi({obj.phys.species},species{i}));
                
                cl_mic = sp_in.CLint_mic;
                fu_mic = sp_in.fuinc_mic;
                cl_hep = sp_in.CLint_hep;
                fu_hep = sp_in.fuinc_hep;
                fub = sp_in.fub;
                bp = sp_in.bp;
                
                lbf = sp_phys.LBF_mLminkg;
                sf = sp_phys.SF_gLiver_kgBW;
                
                if ~isempty(cl_mic)
                    % non-restrictive microsomal clearance
                    clb_mic_nonres = (cl_mic * sf * lbf)./((cl_mic * sf) + lbf);
                    clp_mic_nonres = ((cl_mic * sf * lbf)./((cl_mic * sf) + lbf)) * bp;
                    obj = obj.addPropertyValue('ivive','species',species{i},'type','non-restrictive','CLb_mLminkg',clb_mic_nonres,'CLp_mLminkg',clp_mic_nonres,'system','microsomes');
                    % restrictive microsomal clearance
                    clb_mic_res = ((cl_mic/fu_mic) * sf * fub * lbf)./(((cl_mic/fu_mic) * sf * fub) + lbf);
                    clp_mic_res = (((cl_mic/fu_mic) * sf * fub * lbf)./(((cl_mic/fu_mic) * sf * fub) + lbf)) * bp;
                    obj = obj.addPropertyValue('ivive','species',species{i},'type','restrictive','CLb_mLminkg',clb_mic_res,'CLp_mLminkg',clp_mic_res,'system','microsomes');
                end
                if ~isempty(cl_hep)
                    % non-restrictive hepatocyte clearance
                    clb_hep_nonres = (cl_hep * sf * lbf)./((cl_hep * sf) + lbf);
                    clp_hep_nonres = ((cl_hep * sf * lbf)./((cl_hep * sf) + lbf)) * bp;
                    obj = obj.addPropertyValue('ivive','species',species{i},'type','non-restrictive','CLb_mLminkg',clb_hep_nonres,'CLp_mLminkg',clp_hep_nonres,'system','hepatocytes');
                    % restrictive hepatocyte clearance
                    clb_hep_res = ((cl_hep/fu_hep) * sf * fub * lbf)./(((cl_hep/fu_hep) * sf * fub) + lbf);
                    clp_hep_res = (((cl_hep/fu_hep) * sf * fub * lbf)./(((cl_hep/fu_hep) * sf * fub) + lbf)) * bp;
                    obj = obj.addPropertyValue('ivive','species',species{i},'type','restrictive','CLb_mLminkg',clb_hep_res,'CLp_mLminkg',clp_hep_res,'system','hepatocytes');
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%% calculate results from IVIVE %%%%
            pk_data = obj.pk_kate;
            idx_iv = ismember(strrep(strrep(strrep(lower({pk_data.route}),'-',''),' ',''),'_',''),{'iv' 'ivbolus' 'ivinf' 'ivinfusion'});
            pk_data = pk_data(idx_iv);
            switch calc_opts.iv_data
                case 'discrete'
                    pk_data = pk_data(strcmpi({pk_data.study_design},'discrete'));
                case 'cassette'
                    pk_data = pk_data(strcmpi({pk_data.study_design},'cassette'));
            end
            if ~isempty(pk_data)
                for i = 1:length(species)
                    pk_sp = pk_data(strcmpi({pk_data.species},species{i}));
                    % remove all rows with an empty CL value
                    idx_cl = find(cell2mat(cellfun(@(x) ~isempty(x),{pk_sp.TBdy_CL_mlmin_kg},'UniformOutput',false)));
                    cl = struct('val',[],'cv',[],'n',[],'dose_mgkg',[],'study',{},'strain',{},'matrix',{},'study_design',{});
                    if ~isempty(idx_cl)
                        pk_sp = pk_sp(idx_cl);
                        for j = 1:length(pk_sp)
                            if isempty(pk_sp(j).study), pk_sp(j).study = 'NoStudyNumber'; end
                            if isempty(pk_sp(j).strain), pk_sp(j).strain = 'NoStrain'; end
                            if isempty(pk_sp(j).matrix), pk_sp(j).matrix = 'NoMatrix'; end
                        end
                        study = unique({pk_sp.study});
                        strain = unique({pk_sp.strain});
                        matrix = unique({pk_sp.matrix});
                        count = 0;
                        for xx = 1:length(study)
                           for yy = 1:length(strain)
                              for zz = 1:length(matrix)
                                  pk_spi = pk_sp(strcmpi({pk_sp.study},study{xx})&strcmpi({pk_sp.strain},strain{yy})&strcmpi({pk_sp.matrix},matrix{zz}));
                                  if ~isempty(pk_spi)
                                      [~,idx] = sort([pk_spi.dose_mgkg]);
                                      pk_spi = pk_spi(idx);
                                      binID = zeros(length(pk_spi),1);
                                      for ww = 1:length(pk_spi)
                                          if ww == 1
                                              binID(ww) = 1;
                                          elseif pk_spi(ww).dose_mgkg/pk_spi(ww-1).dose_mgkg >= 1.25
                                              binID(ww) = binID(ww-1)+1;
                                          else
                                              binID(ww) = binID(ww-1);
                                          end
                                      end
                                      bins = unique(binID);
                                      for ww = 1:length(bins)
                                          bin_idx = logical(binID==bins(ww));
                                          count = count + 1;
                                          switch calc_opts.pk_data_op
                                              case 'geometric mean'
                                                  cl(count).val = geomean([pk_spi(bin_idx).TBdy_CL_mlmin_kg]);
                                              case 'arithmetic mean'
                                                  cl(count).val = mean([pk_spi(bin_idx).TBdy_CL_mlmin_kg]);
                                              case 'median'
                                                  cl(count).val = median([pk_spi(bin_idx).TBdy_CL_mlmin_kg]);
                                          end
                                          cl(count).cv = (std([pk_spi(bin_idx).TBdy_CL_mlmin_kg])./mean([pk_spi(bin_idx).TBdy_CL_mlmin_kg])) * 100;
                                          cl(count).n = length([pk_spi(bin_idx).TBdy_CL_mlmin_kg]);
                                          cl(count).dose_mgkg = mean([pk_spi(bin_idx).dose_mgkg]);
                                          cl(count).study = pk_spi(find(bin_idx,1)).study;
                                          cl(count).strain = pk_spi(find(bin_idx,1)).strain;
                                          cl(count).matrix = pk_spi(find(bin_idx,1)).matrix;
                                          cl(count).study_design = pk_spi(find(bin_idx,1)).study_design;
                                      end
                                  end
                              end
                           end
                        end
                    end
                    if ~isempty(cl)
                        ivive_sp = obj.ivive(strcmpi({obj.ivive.species},species{i}));
                        
                        clb_mic_nonres = [];
                        clp_mic_nonres = [];
                        clb_hep_nonres = [];
                        clp_hep_nonres = [];
                        
                        clb_mic_res = [];
                        clp_mic_res = [];
                        clb_hep_res = [];
                        clp_hep_res = [];
                        
                        idx = strcmpi({ivive_sp.type},'non-restrictive') & strcmpi({ivive_sp.system},'microsomes');
                        if any(idx), clb_mic_nonres = ivive_sp(idx).CLb_mLminkg; end
                        
                        idx = strcmpi({ivive_sp.type},'non-restrictive') & strcmpi({ivive_sp.system},'microsomes');
                        if any(idx), clp_mic_nonres = ivive_sp(idx).CLp_mLminkg; end
                        
                        idx = strcmpi({ivive_sp.type},'non-restrictive') & strcmpi({ivive_sp.system},'hepatocytes');
                        if any(idx), clb_hep_nonres = ivive_sp(idx).CLb_mLminkg; end
                        
                        idx = strcmpi({ivive_sp.type},'non-restrictive') & strcmpi({ivive_sp.system},'hepatocytes');
                        if any(idx), clp_hep_nonres = ivive_sp(idx).CLp_mLminkg; end
                        
                        idx = strcmpi({ivive_sp.type},'restrictive') & strcmpi({ivive_sp.system},'microsomes');
                        if any(idx), clb_mic_res = ivive_sp(idx).CLb_mLminkg; end
                        
                        idx = strcmpi({ivive_sp.type},'restrictive') & strcmpi({ivive_sp.system},'microsomes');
                        if any(idx), clp_mic_res = ivive_sp(idx).CLp_mLminkg; end
                        
                        idx = strcmpi({ivive_sp.type},'restrictive') & strcmpi({ivive_sp.system},'hepatocytes');
                        if any(idx), clb_hep_res = ivive_sp(idx).CLb_mLminkg; end
                        
                        idx = strcmpi({ivive_sp.type},'restrictive') & strcmpi({ivive_sp.system},'hepatocytes');
                        if any(idx), clp_hep_res = ivive_sp(idx).CLp_mLminkg; end
                        
                        for zz = 1:length(cl)
                            % only calculate fold error if CL was measured
                            % in plasma or blood
                            switch lower(cl(zz).matrix)
                                case 'blood'
                                    if ~isempty(clb_mic_nonres)
                                        fe_mic_nonres = 10.^(abs(log10(clb_mic_nonres/cl(zz).val))) * sign(log10(clb_mic_nonres/cl(zz).val));
                                        obj = obj.addPropertyValue('ivive_results','species',species{i},'strain',cl(zz).strain,'system','microsomes','type','non-restrictive','CL_mLminkg',clb_mic_nonres,...
                                            'matrix','blood','CLref_mLminkg',cl(zz).val,'fold_error',fe_mic_nonres,'CLref_pctCV',cl(zz).cv,'CLref_n',cl(zz).n,'study',cl(zz).study,'study_design',cl(zz).study_design);
                                    end
                                    
                                    if ~isempty(clb_mic_res)
                                        fe_mic_res = 10.^(abs(log10(clb_mic_res/cl(zz).val))) * sign(log10(clb_mic_res/cl(zz).val));
                                        obj = obj.addPropertyValue('ivive_results','species',species{i},'strain',cl(zz).strain,'system','microsomes','type','restrictive','CL_mLminkg',clb_mic_res,...
                                            'matrix','blood','CLref_mLminkg',cl(zz).val,'fold_error',fe_mic_res,'CLref_pctCV',cl(zz).cv,'CLref_n',cl(zz).n,'study',cl(zz).study,'study_design',cl(zz).study_design);
                                    end
                                    
                                    if ~isempty(clb_hep_nonres)
                                        fe_hep_nonres = 10.^(abs(log10(clb_hep_nonres/cl(zz).val))) * sign(log10(clb_hep_nonres/cl(zz).val));
                                        obj = obj.addPropertyValue('ivive_results','species',species{i},'strain',cl(zz).strain,'system','hepatocytes','type','non-restrictive','CL_mLminkg',clb_hep_nonres,...
                                            'matrix','blood','CLref_mLminkg',cl(zz).val,'fold_error',fe_hep_nonres,'CLref_pctCV',cl(zz).cv,'CLref_n',cl(zz).n,'study',cl(zz).study,'study_design',cl(zz).study_design);
                                    end
                                    
                                    if ~isempty(clb_hep_res)
                                        fe_hep_res = 10.^(abs(log10(clb_hep_res/cl(zz).val))) * sign(log10(clb_hep_res/cl(zz).val));
                                        obj = obj.addPropertyValue('ivive_results','species',species{i},'strain',cl(zz).strain,'system','hepatocytes','type','restrictive','CL_mLminkg',clb_hep_res,...
                                            'matrix','blood','CLref_mLminkg',cl(zz).val,'fold_error',fe_hep_res,'CLref_pctCV',cl(zz).cv,'CLref_n',cl(zz).n,'study',cl(zz).study,'study_design',cl(zz).study_design);
                                    end
                                    
                                case 'plasma'
                                    if ~isempty(clp_mic_nonres)
                                        fe_mic_nonres = 10.^(abs(log10(clp_mic_nonres/cl(zz).val))) * sign(log10(clp_mic_nonres/cl(zz).val));
                                        obj = obj.addPropertyValue('ivive_results','species',species{i},'strain',cl(zz).strain,'system','microsomes','type','non-restrictive','CL_mLminkg',clp_mic_nonres,...
                                            'matrix','plasma','CLref_mLminkg',cl(zz).val,'fold_error',fe_mic_nonres,'CLref_pctCV',cl(zz).cv,'CLref_n',cl(zz).n,'study',cl(zz).study,'study_design',cl(zz).study_design);
                                    end
                                    
                                    if ~isempty(clp_mic_res)
                                        fe_mic_res = 10.^(abs(log10(clp_mic_res/cl(zz).val))) * sign(log10(clp_mic_res/cl(zz).val));
                                        obj = obj.addPropertyValue('ivive_results','species',species{i},'strain',cl(zz).strain,'system','microsomes','type','restrictive','CL_mLminkg',clp_mic_res,...
                                            'matrix','plasma','CLref_mLminkg',cl(zz).val,'fold_error',fe_mic_res,'CLref_pctCV',cl(zz).cv,'CLref_n',cl(zz).n,'study',cl(zz).study,'study_design',cl(zz).study_design);
                                    end
                                    
                                    if ~isempty(clp_hep_nonres)
                                        fe_hep_nonres = 10.^(abs(log10(clp_hep_nonres/cl(zz).val))) * sign(log10(clp_hep_nonres/cl(zz).val));
                                        obj = obj.addPropertyValue('ivive_results','species',species{i},'strain',cl(zz).strain,'system','hepatocytes','type','non-restrictive','CL_mLminkg',clp_hep_nonres,...
                                            'matrix','plasma','CLref_mLminkg',cl(zz).val,'fold_error',fe_hep_nonres,'CLref_pctCV',cl(zz).cv,'CLref_n',cl(zz).n,'study',cl(zz).study,'study_design',cl(zz).study_design);
                                    end
                                    
                                    if ~isempty(clp_hep_res)
                                        fe_hep_res = 10.^(abs(log10(clp_hep_res/cl(zz).val))) * sign(log10(clp_hep_res/cl(zz).val));
                                        obj = obj.addPropertyValue('ivive_results','species',species{i},'strain',cl(zz).strain,'system','hepatocytes','type','restrictive','CL_mLminkg',clp_hep_res,...
                                            'matrix','plasma','CLref_mLminkg',cl(zz).val,'fold_error',fe_hep_res,'CLref_pctCV',cl(zz).cv,'CLref_n',cl(zz).n,'study',cl(zz).study,'study_design',cl(zz).study_design);
                                    end
                            end
                        end
                    end
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
        function obj = LBFscaling(obj,calc_opts)
            pk_data = obj.pk_kate;
            idx_iv = ismember(strrep(strrep(strrep(lower({pk_data.route}),'-',''),' ',''),'_',''),{'iv' 'ivbolus' 'ivinf' 'ivinfusion'});
            pk_data = pk_data(idx_iv);
            switch calc_opts.iv_data
                case 'discrete'
                    pk_data = pk_data(strcmpi({pk_data.study_design},'discrete'));
                case 'cassette'
                    pk_data = pk_data(strcmpi({pk_data.study_design},'cassette'));
            end
            lbf_hu = obj.phys(strcmpi({obj.phys.species},'human')).LBF_mLminkg;
            bp_hu = obj.species_inputs(strcmpi({obj.species_inputs.species},'human')).bp;
            fub_hu = obj.species_inputs(strcmpi({obj.species_inputs.species},'human')).fub;
            for i = 1:length(pk_data)
                switch lower(pk_data(i).matrix)
                    case {'plasma','blood'}
                        sp = pk_data(i).species;
                        lbf_animal = [obj.phys(strcmpi({obj.phys.species},sp)).LBF_mLminkg];
                        bp_animal = [obj.species_inputs(strcmpi({obj.species_inputs.species},sp)).bp];
                        fub_animal = [obj.species_inputs(strcmpi({obj.species_inputs.species},sp)).fub];
                        
                        clb_hu_nonres = [];
                        clp_hu_nonres = [];
                        clb_hu_res = [];
                        clp_hu_res = [];
                        hu_t12_nonres = [];
                        hu_t12_res = [];
                        
                        cl = [];
                        t12 = [];
                        if ~isempty(lbf_animal)
                            if ~isempty(pk_data(i).TBdy_CL_mlmin_kg) && ~isempty(pk_data(i).T1_2_hr)
                                if ~isempty(pk_data(i).TBdy_CL_mlmin_kg)
                                    cl = pk_data(i).TBdy_CL_mlmin_kg;
                                    if strcmpi(pk_data(i).matrix,'blood')
                                        clb_hu_nonres = (cl * lbf_hu)./lbf_animal;
                                        clp_hu_nonres = clb_hu_nonres * bp_hu;

                                        clb_hu_res = ((cl/fub_animal) * lbf_hu * fub_hu)./lbf_animal;
                                        clp_hu_res = clb_hu_res * bp_hu;
                                    else
                                        clb_hu_nonres = ((cl/bp_animal) * lbf_hu)./lbf_animal;
                                        clp_hu_nonres = clb_hu_nonres * bp_hu;

                                        clb_hu_res = ((cl/(fub_animal * bp_animal)) * lbf_hu * fub_hu)./lbf_animal;
                                        clp_hu_res = clb_hu_res * bp_hu;
                                    end
                                end
                                if ~isempty(pk_data(i).T1_2_hr)
                                    t12 = pk_data(i).T1_2_hr;
                                    hu_t12_nonres = (t12 * lbf_animal)./lbf_hu;
                                    hu_t12_res = (t12 * lbf_animal * fub_animal)./(lbf_hu * fub_hu);
                                end
                                obj = obj.addPropertyValue('lbf','cmpd',obj.name,'species',sp,'strain',pk_data(i).strain,'sex',pk_data(i).sex,'matrix',pk_data(i).matrix,...
                                    'CL_mLminkg',cl,'Thalf_hr',t12,'type','non-restrictive','HuCLp_mLminkg',clp_hu_nonres,'HuCLb_mLminkg',clb_hu_nonres,'HuThalf_hr',hu_t12_nonres,...
                                    'animal_id',pk_data(i).animal_id,'study',pk_data(i).study,'study_design',pk_data(i).study_design);

                                obj = obj.addPropertyValue('lbf','cmpd',obj.name,'species',sp,'strain',pk_data(i).strain,'sex',pk_data(i).sex,'matrix',pk_data(i).matrix,...
                                    'CL_mLminkg',cl,'Thalf_hr',t12,'type','restrictive','HuCLp_mLminkg',clp_hu_res,'HuCLb_mLminkg',clb_hu_res,'HuThalf_hr',hu_t12_res,...
                                    'animal_id',pk_data(i).animal_id,'study',pk_data(i).study,'study_design',pk_data(i).study_design);
                            end
                        else
                            warning('Species %s not supported for LBF scaling',sp)
                        end
                end
            end
            sp = unique({obj.lbf.species});
            for i = 1:length(sp)
                
                lbf_nonres = obj.lbf(strcmpi({obj.lbf.species},sp{i}) & strcmpi({obj.lbf.type},'non-restrictive'));
                lbf_res = obj.lbf(strcmpi({obj.lbf.species},sp{i}) & strcmpi({obj.lbf.type},'restrictive'));
                
                clb_hu_nonres = [];
                clp_hu_nonres = [];
                
                clb_hu_res = [];
                clp_hu_res = [];
                
                hu_t12_nonres = [];
                hu_t12_res = [];
                
                switch calc_opts.pk_data_op
                    case 'geometric mean'
                        clb_hu_nonres = geomean([lbf_nonres.HuCLb_mLminkg]);
                        clp_hu_nonres = geomean([lbf_nonres.HuCLp_mLminkg]);
                        hu_t12_nonres = geomean([lbf_nonres.HuThalf_hr]);
                        
                        clb_hu_res = geomean([lbf_res.HuCLb_mLminkg]);
                        clp_hu_res = geomean([lbf_res.HuCLp_mLminkg]);
                        hu_t12_res = geomean([lbf_res.HuThalf_hr]);
                    case 'arithmetic mean'
                        clb_hu_nonres = mean([lbf_nonres.HuCLb_mLminkg]);
                        clp_hu_nonres = mean([lbf_nonres.HuCLp_mLminkg]);
                        hu_t12_nonres = mean([lbf_nonres.HuThalf_hr]);
                        
                        clb_hu_res = mean([lbf_res.HuCLb_mLminkg]);
                        clp_hu_res = mean([lbf_res.HuCLp_mLminkg]);
                        hu_t12_res = mean([lbf_res.HuThalf_hr]);
                    case 'median'
                        clb_hu_nonres = median([lbf_nonres.HuCLb_mLminkg]);
                        clp_hu_nonres = median([lbf_nonres.HuCLp_mLminkg]);
                        hu_t12_nonres = median([lbf_nonres.HuThalf_hr]);
                        
                        clb_hu_res = median([lbf_res.HuCLb_mLminkg]);
                        clp_hu_res = median([lbf_res.HuCLp_mLminkg]);
                        hu_t12_res = median([lbf_res.HuThalf_hr]);
                end
                
                obj = obj.addPropertyValue('lbf','cmpd',obj.name,'species',sp{i},'type','restrictive','HuCLp_mLminkg',clp_hu_res,'HuCLb_mLminkg',clb_hu_res,'HuThalf_hr',hu_t12_res,'study','summary');
                obj = obj.addPropertyValue('lbf','cmpd',obj.name,'species',sp{i},'type','non-restrictive','HuCLp_mLminkg',clp_hu_nonres,'HuCLb_mLminkg',clb_hu_nonres,'HuThalf_hr',hu_t12_nonres,'study','summary');
            end
            % add summary of all species
            if ~isempty(obj.lbf)
                lbf_nonres = obj.lbf(strcmpi({obj.lbf.type},'non-restrictive'));
                lbf_res = obj.lbf(strcmpi({obj.lbf.type},'restrictive'));

                clb_hu_nonres = [];
                clp_hu_nonres = [];

                clb_hu_res = [];
                clp_hu_res = [];

                hu_t12_nonres = [];
                hu_t12_res = [];

                switch calc_opts.pk_data_op
                    case 'geometric mean'
                        clb_hu_nonres = geomean([lbf_nonres.HuCLb_mLminkg]);
                        clp_hu_nonres = geomean([lbf_nonres.HuCLp_mLminkg]);
                        hu_t12_nonres = geomean([lbf_nonres.HuThalf_hr]);

                        clb_hu_res = geomean([lbf_res.HuCLb_mLminkg]);
                        clp_hu_res = geomean([lbf_res.HuCLp_mLminkg]);
                        hu_t12_res = geomean([lbf_res.HuThalf_hr]);
                    case 'arithmetic mean'
                        clb_hu_nonres = mean([lbf_nonres.HuCLb_mLminkg]);
                        clp_hu_nonres = mean([lbf_nonres.HuCLp_mLminkg]);
                        hu_t12_nonres = mean([lbf_nonres.HuThalf_hr]);

                        clb_hu_res = mean([lbf_res.HuCLb_mLminkg]);
                        clp_hu_res = mean([lbf_res.HuCLp_mLminkg]);
                        hu_t12_res = mean([lbf_res.HuThalf_hr]);
                    case 'median'
                        clb_hu_nonres = median([lbf_nonres.HuCLb_mLminkg]);
                        clp_hu_nonres = median([lbf_nonres.HuCLp_mLminkg]);
                        hu_t12_nonres = median([lbf_nonres.HuThalf_hr]);

                        clb_hu_res = median([lbf_res.HuCLb_mLminkg]);
                        clp_hu_res = median([lbf_res.HuCLp_mLminkg]);
                        hu_t12_res = median([lbf_res.HuThalf_hr]);
                end

                obj = obj.addPropertyValue('lbf','cmpd',obj.name,'species','all','type','restrictive','HuCLp_mLminkg',clp_hu_res,'HuCLb_mLminkg',clb_hu_res,'HuThalf_hr',hu_t12_res,...
                    'study','summary','study_design',calc_opts.pk_data_op);
                obj = obj.addPropertyValue('lbf','cmpd',obj.name,'species','all','type','non-restrictive','HuCLp_mLminkg',clp_hu_nonres,'HuCLb_mLminkg',clb_hu_nonres,'HuThalf_hr',hu_t12_nonres,...
                    'study','summary','study_design',calc_opts.pk_data_op);
            end
            
        end
        
        function obj = SingleSpeciesAllometry(obj,calc_opts)
            pk_data = obj.pk_kate;
            idx_iv = ismember(strrep(strrep(strrep(lower({pk_data.route}),'-',''),' ',''),'_',''),{'iv' 'ivbolus' 'ivinf' 'ivinfusion'});
            pk_data = pk_data(idx_iv);
            idx_matrix = ismember(lower({pk_data.matrix}),{'plasma' 'blood'});
            pk_data = pk_data(idx_matrix);
            switch calc_opts.iv_data
                case 'discrete'
                    pk_data = pk_data(strcmpi({pk_data.study_design},'discrete'));
                case 'cassette'
                    pk_data = pk_data(strcmpi({pk_data.study_design},'cassette'));
            end
            bw_hu = obj.phys(strcmpi({obj.phys.species},'human')).BW_kg;
            bp_hu = obj.species_inputs(strcmpi({obj.species_inputs.species},'human')).bp;
            fub_hu = obj.species_inputs(strcmpi({obj.species_inputs.species},'human')).fub;
            for i = 1:length(pk_data)
                switch lower(pk_data(i).matrix)
                    case {'plasma','blood'}
                        sp = pk_data(i).species;
                        bw_animal = [obj.phys(strcmpi({obj.phys.species},sp)).BW_kg];
                        bp_animal = [obj.species_inputs(strcmpi({obj.species_inputs.species},sp)).bp];
                        fub_animal = [obj.species_inputs(strcmpi({obj.species_inputs.species},sp)).fub];
                        
                        clb_hu_nonres = [];
                        clp_hu_nonres = [];
                        clb_hu_res = [];
                        clp_hu_res = [];
                        vssb_hu_nonres = [];
                        vssp_hu_nonres = [];
                        vssb_hu_res = [];
                        vssp_hu_res = [];
                        hu_t12_nonres = [];
                        hu_t12_res = [];
                        
                        clb = [];
                        clp = [];
                        vssb = [];
                        vssp = [];
                        
                        if ~isempty(bw_animal)
                            if ~isempty(pk_data(i).TBdy_CL_mlmin_kg) && ~isempty(pk_data(i).T1_2_hr)
                                if ~isempty(pk_data(i).TBdy_CL_mlmin_kg)
                                    cl = pk_data(i).TBdy_CL_mlmin_kg * bw_animal;
                                    if strcmpi(pk_data(i).matrix,'blood')
                                        clb_hu_nonres = (cl * ((bw_hu/bw_animal)^0.75)) / bw_hu;
                                        clp_hu_nonres = clb_hu_nonres * bp_hu;

                                        clb_hu_res = ((cl/fub_animal) * ((bw_hu/bw_animal)^0.75) * fub_hu) / bw_hu;
                                        clp_hu_res = clb_hu_res * bp_hu;
                                        
                                        clb = cl;
                                        clp = cl * bp_animal;
                                    else
                                        clb_hu_nonres = ((cl/bp_animal) * ((bw_hu/bw_animal)^0.75)) / bw_hu;
                                        clp_hu_nonres = clb_hu_nonres * bp_hu;

                                        clb_hu_res = ((cl/(fub_animal*bp_animal)) * ((bw_hu/bw_animal)^0.75) * fub_hu) / bw_hu;
                                        clp_hu_res = clb_hu_res * bp_hu;
                                        
                                        clb = cl / bp_animal;
                                        clp = cl;
                                    end
                                end
                                if ~isempty(pk_data(i).Vdss_Lkg)
                                    vss = pk_data(i).Vdss_Lkg;
                                    if strcmpi(pk_data(i).matrix,'blood')
                                        vssb_hu_nonres = vss;
                                        vssp_hu_nonres = vssb_hu_nonres * bp_hu;
                                        
                                        vssb_hu_res = (vss/fub_animal) * fub_hu;
                                        vssp_hu_res = vssb_hu_res * bp_hu;
                                        
                                        vssb = vss;
                                        vssp = vss * bp_animal;
                                    else
                                        vssb_hu_nonres = vss/bp_animal;
                                        vssp_hu_nonres = vssb_hu_nonres * bp_hu;
                                        
                                        vssb_hu_res = (vss/(bp_animal*fub_animal)) * fub_hu;
                                        vssp_hu_res = vssb_hu_res * bp_hu;
                                        
                                        vssb = vss / bp_animal;
                                        vssp = vss;
                                    end
                                end
                                if ~isempty(vssb_hu_nonres) && ~isempty(clb_hu_nonres)
                                    hu_t12_nonres = (vssb_hu_nonres*log(2))/clb_hu_nonres;
                                end
                                if ~isempty(vssb_hu_res) && ~isempty(clb_hu_res)
                                    hu_t12_res = (vssb_hu_res*log(2))/clb_hu_res;
                                end
                                
                                obj = obj.addPropertyValue('single_allometry','cmpd',obj.name,'species',sp,'strain',pk_data(i).strain,'sex',pk_data(i).sex,'matrix',pk_data(i).matrix,...
                                    'CLb_mLminkg',clb/bw_animal,'CLp_mLminkg',clp/bw_animal,'Vss_b_Lkg',vssb,'Vss_p_Lkg',vssp,'type','non-restrictive',...
                                    'HuCLp_mLminkg',clp_hu_nonres,'HuCLb_mLminkg',clb_hu_nonres,'HuVss_b_Lkg',vssb_hu_nonres,'HuVss_p_Lkg',vssp_hu_nonres,'HuThalf_hr',hu_t12_nonres,...
                                    'animal_id',pk_data(i).animal_id,'study',pk_data(i).study,'study_design',pk_data(i).study_design);

                                obj = obj.addPropertyValue('single_allometry','cmpd',obj.name,'species',sp,'strain',pk_data(i).strain,'sex',pk_data(i).sex,'matrix',pk_data(i).matrix,...
                                    'CLb_mLminkg',clb/bw_animal,'CLp_mLminkg',clp/bw_animal,'Vss_b_Lkg',vssb,'Vss_p_Lkg',vssp,'type','restrictive',...
                                    'HuCLp_mLminkg',clp_hu_res,'HuCLb_mLminkg',clb_hu_res,'HuVss_b_Lkg',vssb_hu_res,'HuVss_p_Lkg',vssp_hu_res,'HuThalf_hr',hu_t12_res,...
                                    'animal_id',pk_data(i).animal_id,'study',pk_data(i).study,'study_design',pk_data(i).study_design);
                            end
                        else
                            warning('Species %s not supported for allometric scaling',sp)
                        end
                end
            end
            
            allom = obj.single_allometry;
            if ~isempty(allom)
                sp = unique({allom.species});
                for i = 1:length(sp)
                    allom_i = allom(strcmpi({allom.species},sp{i}));
                    
                    clb = [allom_i(~isnan([allom_i.CLb_mLminkg]) & strcmpi({allom_i.type},'non-restrictive')).CLb_mLminkg];
                    clp = [allom_i(~isnan([allom_i.CLp_mLminkg]) & strcmpi({allom_i.type},'non-restrictive')).CLp_mLminkg];
                    vssb = [allom_i(~isnan([allom_i.Vss_b_Lkg]) & strcmpi({allom_i.type},'non-restrictive')).Vss_b_Lkg];
                    vssp = [allom_i(~isnan([allom_i.Vss_p_Lkg]) & strcmpi({allom_i.type},'non-restrictive')).Vss_p_Lkg];
                    
                    hu_clb_nonres = [allom_i(~isnan([allom_i.HuCLb_mLminkg]) & strcmpi({allom_i.type},'non-restrictive')).HuCLb_mLminkg];
                    hu_clp_nonres = [allom_i(~isnan([allom_i.HuCLp_mLminkg]) & strcmpi({allom_i.type},'non-restrictive')).HuCLp_mLminkg];
                    hu_vssb_nonres = [allom_i(~isnan([allom_i.HuVss_b_Lkg]) & strcmpi({allom_i.type},'non-restrictive')).HuVss_b_Lkg];
                    hu_vssp_nonres = [allom_i(~isnan([allom_i.HuVss_p_Lkg]) & strcmpi({allom_i.type},'non-restrictive')).HuVss_p_Lkg];
                    hu_t12_nonres = [allom_i(~isnan([allom_i.HuThalf_hr]) & strcmpi({allom_i.type},'non-restrictive')).HuThalf_hr];
                    
                    hu_clb_res = [allom_i(~isnan([allom_i.HuCLb_mLminkg]) & strcmpi({allom_i.type},'restrictive')).HuCLb_mLminkg];
                    hu_clp_res = [allom_i(~isnan([allom_i.HuCLp_mLminkg]) & strcmpi({allom_i.type},'restrictive')).HuCLp_mLminkg];
                    hu_vssb_res = [allom_i(~isnan([allom_i.HuVss_b_Lkg]) & strcmpi({allom_i.type},'restrictive')).HuVss_b_Lkg];
                    hu_vssp_res = [allom_i(~isnan([allom_i.HuVss_p_Lkg]) & strcmpi({allom_i.type},'restrictive')).HuVss_p_Lkg];
                    hu_t12_res = [allom_i(~isnan([allom_i.HuThalf_hr]) & strcmpi({allom_i.type},'restrictive')).HuThalf_hr];
                    
                    switch lower(calc_opts.pk_data_op)
                        case 'geometric mean'
                            clb = geomean(clb);
                            clp = geomean(clp);
                            vssb = geomean(vssb);
                            vssp = geomean(vssp);
                            
                            hu_clb_nonres = geomean(hu_clb_nonres);
                            hu_clp_nonres = geomean(hu_clp_nonres);
                            hu_vssb_nonres = geomean(hu_vssb_nonres);
                            hu_vssp_nonres = geomean(hu_vssp_nonres);
                            hu_t12_nonres = geomean(hu_t12_nonres);
                            
                            hu_clb_res = geomean(hu_clb_res);
                            hu_clp_res = geomean(hu_clp_res);
                            hu_vssb_res = geomean(hu_vssb_res);
                            hu_vssp_res = geomean(hu_vssp_res);
                            hu_t12_res = geomean(hu_t12_res);
                        case 'arithmetic mean'
                            clb = mean(clb);
                            clp = mean(clp);
                            vssb = mean(vssb);
                            vssp = mean(vssp);
                            
                            hu_clb_nonres = mean(hu_clb_nonres);
                            hu_clp_nonres = mean(hu_clp_nonres);
                            hu_vssb_nonres = mean(hu_vssb_nonres);
                            hu_vssp_nonres = mean(hu_vssp_nonres);
                            hu_t12_nonres = mean(hu_t12_nonres);
                            
                            hu_clb_res = mean(hu_clb_res);
                            hu_clp_res = mean(hu_clp_res);
                            hu_vssb_res = mean(hu_vssb_res);
                            hu_vssp_res = mean(hu_vssp_res);
                            hu_t12_res = mean(hu_t12_res);
                        case 'median'
                            clb = median(clb);
                            clp = median(clp);
                            vssb = median(vssb);
                            vssp = median(vssp);
                            
                            hu_clb_nonres = median(hu_clb_nonres);
                            hu_clp_nonres = median(hu_clp_nonres);
                            hu_vssb_nonres = median(hu_vssb_nonres);
                            hu_vssp_nonres = median(hu_vssp_nonres);
                            hu_t12_nonres = median(hu_t12_nonres);
                            
                            hu_clb_res = median(hu_clb_res);
                            hu_clp_res = median(hu_clp_res);
                            hu_vssb_res = median(hu_vssb_res);
                            hu_vssp_res = median(hu_vssp_res);
                            hu_t12_res = median(hu_t12_res);
                    end
                    obj = obj.addPropertyValue('single_allometry','cmpd',obj.name,'species',sp{i},...
                        'CLb_mLminkg',clb,'CLp_mLminkg',clp,'Vss_b_Lkg',vssb,'Vss_p_Lkg',vssp,'type','restrictive',...
                        'HuCLp_mLminkg',hu_clp_res,'HuCLb_mLminkg',hu_clb_res,'HuVss_b_Lkg',hu_vssb_res,'HuVss_p_Lkg',hu_vssp_res,'HuThalf_hr',hu_t12_res,...
                        'study','summary','study_design',calc_opts.pk_data_op);
                    
                    obj = obj.addPropertyValue('single_allometry','cmpd',obj.name,'species',sp{i},...
                        'CLb_mLminkg',clb,'CLp_mLminkg',clp,'Vss_b_Lkg',vssb,'Vss_p_Lkg',vssp,'type','non-restrictive',...
                        'HuCLp_mLminkg',hu_clp_nonres,'HuCLb_mLminkg',hu_clb_nonres,'HuVss_b_Lkg',hu_vssb_nonres,'HuVss_p_Lkg',hu_vssp_nonres,'HuThalf_hr',hu_t12_nonres,...
                        'study','summary','study_design',calc_opts.pk_data_op);
                end
                
                allom_i = allom(~strcmpi({allom.study},'summary'));
                
                clb = [allom_i(~isnan([allom_i.CLb_mLminkg]) & strcmpi({allom_i.type},'non-restrictive')).CLb_mLminkg];
                clp = [allom_i(~isnan([allom_i.CLp_mLminkg]) & strcmpi({allom_i.type},'non-restrictive')).CLp_mLminkg];
                vssb = [allom_i(~isnan([allom_i.Vss_b_Lkg]) & strcmpi({allom_i.type},'non-restrictive')).Vss_b_Lkg];
                vssp = [allom_i(~isnan([allom_i.Vss_p_Lkg]) & strcmpi({allom_i.type},'non-restrictive')).Vss_p_Lkg];

                hu_clb_nonres = [allom_i(~isnan([allom_i.HuCLb_mLminkg]) & strcmpi({allom_i.type},'non-restrictive')).HuCLb_mLminkg];
                hu_clp_nonres = [allom_i(~isnan([allom_i.HuCLp_mLminkg]) & strcmpi({allom_i.type},'non-restrictive')).HuCLp_mLminkg];
                hu_vssb_nonres = [allom_i(~isnan([allom_i.HuVss_b_Lkg]) & strcmpi({allom_i.type},'non-restrictive')).HuVss_b_Lkg];
                hu_vssp_nonres = [allom_i(~isnan([allom_i.HuVss_p_Lkg]) & strcmpi({allom_i.type},'non-restrictive')).HuVss_p_Lkg];
                hu_t12_nonres = [allom_i(~isnan([allom_i.HuThalf_hr]) & strcmpi({allom_i.type},'non-restrictive')).HuThalf_hr];

                hu_clb_res = [allom_i(~isnan([allom_i.HuCLb_mLminkg]) & strcmpi({allom_i.type},'restrictive')).HuCLb_mLminkg];
                hu_clp_res = [allom_i(~isnan([allom_i.HuCLp_mLminkg]) & strcmpi({allom_i.type},'restrictive')).HuCLp_mLminkg];
                hu_vssb_res = [allom_i(~isnan([allom_i.HuVss_b_Lkg]) & strcmpi({allom_i.type},'restrictive')).HuVss_b_Lkg];
                hu_vssp_res = [allom_i(~isnan([allom_i.HuVss_p_Lkg]) & strcmpi({allom_i.type},'restrictive')).HuVss_p_Lkg];
                hu_t12_res = [allom_i(~isnan([allom_i.HuThalf_hr]) & strcmpi({allom_i.type},'restrictive')).HuThalf_hr];
                
                switch lower(calc_opts.pk_data_op)
                    case 'geometric mean'
                        clb = geomean(clb);
                        clp = geomean(clp);
                        vssb = geomean(vssb);
                        vssp = geomean(vssp);

                        hu_clb_nonres = geomean(hu_clb_nonres);
                        hu_clp_nonres = geomean(hu_clp_nonres);
                        hu_vssb_nonres = geomean(hu_vssb_nonres);
                        hu_vssp_nonres = geomean(hu_vssp_nonres);
                        hu_t12_nonres = geomean(hu_t12_nonres);

                        hu_clb_res = geomean(hu_clb_res);
                        hu_clp_res = geomean(hu_clp_res);
                        hu_vssb_res = geomean(hu_vssb_res);
                        hu_vssp_res = geomean(hu_vssp_res);
                        hu_t12_res = geomean(hu_t12_res);
                    case 'arithmetic mean'
                        clb = mean(clb);
                        clp = mean(clp);
                        vssb = mean(vssb);
                        vssp = mean(vssp);

                        hu_clb_nonres = mean(hu_clb_nonres);
                        hu_clp_nonres = mean(hu_clp_nonres);
                        hu_vssb_nonres = mean(hu_vssb_nonres);
                        hu_vssp_nonres = mean(hu_vssp_nonres);
                        hu_t12_nonres = mean(hu_t12_nonres);

                        hu_clb_res = mean(hu_clb_res);
                        hu_clp_res = mean(hu_clp_res);
                        hu_vssb_res = mean(hu_vssb_res);
                        hu_vssp_res = mean(hu_vssp_res);
                        hu_t12_res = mean(hu_t12_res);
                    case 'median'
                        clb = median(clb);
                        clp = median(clp);
                        vssb = median(vssb);
                        vssp = median(vssp);

                        hu_clb_nonres = median(hu_clb_nonres);
                        hu_clp_nonres = median(hu_clp_nonres);
                        hu_vssb_nonres = median(hu_vssb_nonres);
                        hu_vssp_nonres = median(hu_vssp_nonres);
                        hu_t12_nonres = median(hu_t12_nonres);

                        hu_clb_res = median(hu_clb_res);
                        hu_clp_res = median(hu_clp_res);
                        hu_vssb_res = median(hu_vssb_res);
                        hu_vssp_res = median(hu_vssp_res);
                        hu_t12_res = median(hu_t12_res);
                end
                
                obj = obj.addPropertyValue('single_allometry','cmpd',obj.name,'species','all',...
                    'CLb_mLminkg',clb,'CLp_mLminkg',clp,'Vss_b_Lkg',vssb,'Vss_p_Lkg',vssp,'type','restrictive',...
                    'HuCLp_mLminkg',hu_clp_res,'HuCLb_mLminkg',hu_clb_res,'HuVss_b_Lkg',hu_vssb_res,'HuVss_p_Lkg',hu_vssp_res,'HuThalf_hr',hu_t12_res,...
                    'study','summary','study_design',calc_opts.pk_data_op);

                obj = obj.addPropertyValue('single_allometry','cmpd',obj.name,'species','all',...
                    'CLb_mLminkg',clb,'CLp_mLminkg',clp,'Vss_b_Lkg',vssb,'Vss_p_Lkg',vssp,'type','non-restrictive',...
                    'HuCLp_mLminkg',hu_clp_nonres,'HuCLb_mLminkg',hu_clb_nonres,'HuVss_b_Lkg',hu_vssb_nonres,'HuVss_p_Lkg',hu_vssp_nonres,'HuThalf_hr',hu_t12_nonres,...
                    'study','summary','study_design',calc_opts.pk_data_op);
                
            end
        end
        
        function obj = BuildAllometryData(obj,calc_opts)
            pk_data = obj.pk_kate;
            idx_iv = ismember(strrep(strrep(strrep(lower({pk_data.route}),'-',''),' ',''),'_',''),{'iv' 'ivbolus' 'ivinf' 'ivinfusion'});
            pk_data = pk_data(idx_iv);
            switch calc_opts.iv_data
                case 'discrete'
                    pk_data = pk_data(strcmpi({pk_data.study_design},'discrete'));
                case 'cassette'
                    pk_data = pk_data(strcmpi({pk_data.study_design},'cassette'));
            end
            for i = 1:length(pk_data)
                switch lower(pk_data(i).matrix)
                    case {'blood','plasma'}
                        sp = pk_data(i).species;
                        bw_animal = [obj.phys(strcmpi({obj.phys.species},sp)).BW_kg];
                        bp_animal = [obj.species_inputs(strcmpi({obj.species_inputs.species},sp)).bp];
                        fub_animal = [obj.species_inputs(strcmpi({obj.species_inputs.species},sp)).fub];
                        fup_animal = [obj.species_inputs(strcmpi({obj.species_inputs.species},sp)).fup];
                        brw_animal = [obj.phys(strcmpi({obj.phys.species},sp)).BrW_g];
                        mlp_animal = [obj.phys(strcmpi({obj.phys.species},sp)).MLP_yrs];
                        
                        cl = [];
                        vss = [];
                        if ~isempty(bw_animal)
                            if strcmpi(pk_data(i).matrix,'blood')
                                if ~isempty(pk_data(i).TBdy_CL_mlmin_kg)
                                    cl = pk_data(i).TBdy_CL_mlmin_kg;
                                end
                                if ~isempty(pk_data(i).Vdss_Lkg)
                                    vss = pk_data(i).Vdss_Lkg;
                                end
                                obj = obj.addPropertyValue('allometry_data','cmpd',obj.name,'species',sp,'BW_kg',bw_animal,'BrW_g',brw_animal,'MLP_yrs',mlp_animal,'matrix',pk_data(i).matrix,...
                                    'CL_mLminkg',cl,'Vss_Lkg',vss,'CLb_mLmin',cl * bw_animal,'CLp_mLmin',cl * bw_animal * bp_animal,'CLbU_mLmin',(cl * bw_animal) / fub_animal,'CLpU_mLmin',(cl * bw_animal * bp_animal) / fup_animal,...
                                    'Vss_b_L',vss * bw_animal,'Vss_p_L',vss * bw_animal * bp_animal,'VssU_b_L',vss * bw_animal * fub_animal,'VssU_p_L',vss * bw_animal * bp_animal * fup_animal,...
                                    'CLb_x_MLP',cl * bw_animal * mlp_animal,'CLp_x_MLP',cl * bw_animal * bp_animal * mlp_animal,'CLb_x_BrW',cl * bw_animal * brw_animal,'CLp_x_BrW',cl * bw_animal * bp_animal * brw_animal,...
                                    'CLbU_x_MLP',(cl * bw_animal * mlp_animal) / fub_animal,'CLpU_x_MLP',(cl * bw_animal * bp_animal * mlp_animal) / fup_animal,'CLbU_x_BrW',(cl * bw_animal * brw_animal) / fub_animal,'CLpU_x_BrW',(cl * bw_animal * bp_animal * brw_animal) / fup_animal);
                            else
                                if ~isempty(pk_data(i).TBdy_CL_mlmin_kg)
                                    cl = pk_data(i).TBdy_CL_mlmin_kg;
                                else
                                    cl = NaN;
                                end
                                if ~isempty(pk_data(i).Vdss_Lkg)
                                    vss = pk_data(i).Vdss_Lkg;
                                else
                                    vss = NaN;
                                end
                                if ~isnan(cl) && ~isnan(vss)
                                    obj = obj.addPropertyValue('allometry_data','cmpd',obj.name,'species',sp,'BW_kg',bw_animal,'BrW_g',brw_animal,'MLP_yrs',mlp_animal,'matrix',pk_data(i).matrix,...
                                        'CL_mLminkg',cl,'Vss_Lkg',vss,'CLb_mLmin',(cl/bp_animal) * bw_animal,'CLp_mLmin',cl * bw_animal,'CLbU_mLmin',((cl/bp_animal) * bw_animal) / fub_animal,'CLpU_mLmin',(cl * bw_animal) / fup_animal,...
                                        'Vss_b_L',(vss/bp_animal) * bw_animal,'Vss_p_L',vss * bw_animal,'VssU_b_L',((vss/bp_animal) * bw_animal) / fub_animal,'VssU_p_L',(vss * bw_animal) / fup_animal,...
                                        'CLb_x_MLP',(cl/bp_animal) * bw_animal * mlp_animal,'CLp_x_MLP',cl * bw_animal * mlp_animal,'CLb_x_BrW',(cl/bp_animal) * bw_animal * brw_animal,'CLp_x_BrW',cl * bw_animal * brw_animal,...
                                        'CLbU_x_MLP',((cl/bp_animal) * bw_animal * mlp_animal) / fub_animal,'CLpU_x_MLP',(cl * bw_animal * mlp_animal) / fup_animal,'CLbU_x_BrW',((cl/bp_animal) * bw_animal * brw_animal) / fub_animal,'CLpU_x_BrW',(cl * bw_animal * brw_animal) / fup_animal);
                                end
                            end
                        end
                end
            end
        end
        
        function obj = MultiSpeciesAllometry(obj)
            allom = obj.allometry_data;
            idx = arrayfun(@(allom) ~isempty(allom.CLbU_mLmin),allom);
            allom = allom(idx);
            sp = unique({allom.species});
            
            run = 0;
            % run multi-species allometry if 2 species are available, ONLY
            % IF dog or minipig is one of the species
            if length(sp) == 2
                if any(ismember(strrep(lower(sp),' ',''),{'minipig','pig','dog'}))
                    run = 1;
                end
            elseif length(sp) >= 3 % run if 3 or more species are available
                run = 1;
            end
            
            if run  % initially allow allometry calculations if there are two or more species
                
                bw_human = obj.phys(strcmpi({obj.phys.species},'human')).BW_kg;
                brw_human = obj.phys(strcmpi({obj.phys.species},'human')).BrW_g;
                mlp_human = obj.phys(strcmpi({obj.phys.species},'human')).MLP_yrs;
                
                fub_human = obj.species_inputs(strcmpi({obj.species_inputs.species},'human')).fub;
                bp_human = obj.species_inputs(strcmpi({obj.species_inputs.species},'human')).bp;
                
                log_bw = log10([allom.BW_kg])';
                
                log_CL = log10([allom.CLb_mLmin])';
                
                log_CLu = log10([allom.CLbU_mLmin])';
                
                log_CL_MLP = log10([allom.CLb_x_MLP])';
                
                log_CLu_MLP = log10([allom.CLbU_x_MLP])';
                
                log_CL_BrW = log10([allom.CLb_x_BrW])';
                
                log_CLu_BrW = log10([allom.CLbU_x_BrW])';
                
                log_Vss = log10([allom.Vss_b_L])';
                
                log_VssU = log10([allom.VssU_b_L])';
                
                p_cl_simple = polyfit(log_bw,log_CL,1);
                pred_cl_simple = 10.^polyval(p_cl_simple,log10(bw_human)) / bw_human;
                cor = corrcoef(log_CL,polyval(p_cl_simple,log_bw));
                r2_cl_simple = cor(2,1).^2;
                
                obj = obj.addPropertyValue('multi_allometry','cmpd',obj.name,'parameter','CL_mLminkg','method','Simple','type','non-restrictive',...
                    'value_blood',pred_cl_simple,'value_plasma',pred_cl_simple * bp_human,'R2_fit',r2_cl_simple,'a',10^p_cl_simple(2),'b',p_cl_simple(1),...
                    'equation',sprintf('CL = %.*f x BW^%.*f',4,10^p_cl_simple(2),4,p_cl_simple(1)),'preferred_cl_res',0,'preferred_cl_nonres',0);
                
                p_clu_simple = polyfit(log_bw,log_CLu,1);
                pred_clu_simple = (10.^polyval(p_clu_simple,log10(bw_human)) / bw_human) * fub_human;
                cor = corrcoef(log_CLu,polyval(p_clu_simple,log_bw));
                r2_clu_simple = cor(2,1).^2;
                
                obj = obj.addPropertyValue('multi_allometry','cmpd',obj.name,'parameter','CL_mLminkg','method','Simple','type','restrictive',...
                    'value_blood',pred_clu_simple,'value_plasma',pred_clu_simple * bp_human,'R2_fit',r2_clu_simple,'a',10^p_clu_simple(2),'b',p_clu_simple(1),...
                    'equation',sprintf('CL = %.*f x BW^%.*f',4,10^p_clu_simple(2),4,p_clu_simple(1)),'preferred_cl_res',0,'preferred_cl_nonres',0);
                
                p_cl_mlp = polyfit(log_bw,log_CL_MLP,1);
                pred_cl_mlp = 10.^polyval(p_cl_mlp,log10(bw_human)) / mlp_human / bw_human;
                cor = corrcoef(log_CL_MLP,polyval(p_cl_mlp,log_bw));
                r2_cl_mlp = cor(2,1).^2;
                
                obj = obj.addPropertyValue('multi_allometry','cmpd',obj.name,'parameter','CL_mLminkg','method','MLP','type','non-restrictive',...
                    'value_blood',pred_cl_mlp,'value_plasma',pred_cl_mlp * bp_human,'R2_fit',r2_cl_mlp,'a',10^p_cl_mlp(2),'b',p_cl_mlp(1),...
                    'equation',sprintf('CL x MLP = %.*f x BW^%.*f',4,10^p_cl_mlp(2),4,p_cl_mlp(1)),'preferred_cl_res',0,'preferred_cl_nonres',0);
                
                p_clu_mlp = polyfit(log_bw,log_CLu_MLP,1);
                pred_clu_mlp = (10.^polyval(p_clu_mlp,log10(bw_human)) / mlp_human / bw_human) * fub_human;
                cor = corrcoef(log_CLu_MLP,polyval(p_clu_mlp,log_bw));
                r2_clu_mlp = cor(2,1).^2;
                
                obj = obj.addPropertyValue('multi_allometry','cmpd',obj.name,'parameter','CL_mLminkg','method','MLP','type','restrictive',...
                    'value_blood',pred_clu_mlp,'value_plasma',pred_clu_mlp * bp_human,'R2_fit',r2_clu_mlp,'a',10^p_clu_mlp(2),'b',p_clu_mlp(1),...
                    'equation',sprintf('CL x MLP = %.*f x BW^%.*f',4,10^p_clu_mlp(2),4,p_clu_mlp(1)),'preferred_cl_res',0,'preferred_cl_nonres',0);
                
                p_cl_brw = polyfit(log_bw,log_CL_BrW,1);
                pred_cl_brw = 10.^polyval(p_cl_brw,log10(bw_human)) / brw_human / bw_human;
                cor = corrcoef(log_CL_BrW,polyval(p_cl_brw,log_bw));
                r2_cl_brw = cor(2,1).^2;
                
                obj = obj.addPropertyValue('multi_allometry','cmpd',obj.name,'parameter','CL_mLminkg','method','BrW','type','non-restrictive',...
                    'value_blood',pred_cl_brw,'value_plasma',pred_cl_brw * bp_human,'R2_fit',r2_cl_brw,'a',10^p_cl_brw(2),'b',p_cl_brw(1),...
                    'equation',sprintf('CL x BrW = %.*f x BW^%.*f',4,10^p_cl_brw(2),4,p_cl_brw(1)),'preferred_cl_res',0,'preferred_cl_nonres',0);
                
                p_clu_brw = polyfit(log_bw,log_CLu_BrW,1);
                pred_clu_brw = (10.^polyval(p_clu_brw,log10(bw_human)) / brw_human / bw_human) * fub_human;
                cor = corrcoef(log_CLu_BrW,polyval(p_clu_brw,log_bw));
                r2_clu_brw = cor(2,1).^2;
                
                obj = obj.addPropertyValue('multi_allometry','cmpd',obj.name,'parameter','CL_mLminkg','method','BrW','type','restrictive',...
                    'value_blood',pred_clu_brw,'value_plasma',pred_clu_brw * bp_human,'R2_fit',r2_clu_brw,'a',10^p_clu_brw(2),'b',p_clu_brw(1),...
                    'equation',sprintf('CL x BrW = %.*f x BW^%.*f',4,10^p_clu_brw(2),4,p_clu_brw(1)),'preferred_cl_res',0,'preferred_cl_nonres',0);
                
                p_vss = polyfit(log_bw,log_Vss,1);
                pred_vss = 10.^polyval(p_vss,log10(bw_human)) / bw_human;
                cor = corrcoef(log_Vss,polyval(p_vss,log_bw));
                r2_vss = cor(2,1).^2;
                
                obj = obj.addPropertyValue('multi_allometry','cmpd',obj.name,'parameter','Vss_Lkg','method','Simple','type','non-restrictive',...
                    'value_blood',pred_vss,'value_plasma',pred_vss * bp_human,'R2_fit',r2_vss,'a',10^p_vss(2),'b',p_vss(1),...
                    'equation',sprintf('Vss = %.*f x BW^%.*f',4,10^p_vss(2),4,p_vss(1)),'preferred_cl_res',0,'preferred_cl_nonres',0);
                
                p_vssu = polyfit(log_bw,log_VssU,1);
                pred_vssu = (10.^polyval(p_vss,log10(bw_human)) / bw_human) * fub_human;
                cor = corrcoef(log_VssU,polyval(p_vssu,log_bw));
                r2_vssu = cor(2,1).^2;
                
                obj = obj.addPropertyValue('multi_allometry','cmpd',obj.name,'parameter','Vss_Lkg','method','Simple','type','restrictive',...
                    'value_blood',pred_vssu,'value_plasma',pred_vssu * bp_human,'R2_fit',r2_vssu,'a',10^p_vssu(2),'b',p_vssu(1),...
                    'equation',sprintf('Vss = %.*f x BW^%.*f',4,10^p_vssu(2),4,p_vssu(1)),'preferred_cl_res',0,'preferred_cl_nonres',0);
                
                res = obj.multi_allometry(strcmpi({obj.multi_allometry.type},'restrictive') & strcmpi({obj.multi_allometry.parameter},'CL_mLminkg'));
                res_exp = res(strcmpi({res.method},'Simple')).b;
                
                nonres = obj.multi_allometry(strcmpi({obj.multi_allometry.type},'non-restrictive') & strcmpi({obj.multi_allometry.parameter},'CL_mLminkg'));
                nonres_exp = nonres(strcmpi({nonres.method},'Simple')).b;
                
                % use Rule of Exponents (ROE) to define preferred method.
                % ROE reference is in Mahmood, I and J.D. Balian (1996).
                % Xenobiotica 26(9): 887-895. Difference here is the
                % authors put a lower limit on simple allometry of 0.55 -
                % in our case, we assume anything below 0.70 should be
                % simple allometry
                if round(res_exp,2) <= 0.70
                    idx_pref = strcmpi({obj.multi_allometry.type},'restrictive') & strcmpi({obj.multi_allometry.parameter},'CL_mLminkg') & strcmpi({obj.multi_allometry.method},'Simple');
                    obj.multi_allometry(idx_pref).preferred_cl_res = 1;
                elseif round(res_exp,2) > 0.70 && round(res_exp,2) <=1.00
                    idx_pref = strcmpi({obj.multi_allometry.type},'restrictive') & strcmpi({obj.multi_allometry.parameter},'CL_mLminkg') & strcmpi({obj.multi_allometry.method},'MLP');
                    obj.multi_allometry(idx_pref).preferred_cl_res = 1;
                elseif round(res_exp,2) > 1.00
                    idx_pref = strcmpi({obj.multi_allometry.type},'restrictive') & strcmpi({obj.multi_allometry.parameter},'CL_mLminkg') & strcmpi({obj.multi_allometry.method},'BrW');
                    obj.multi_allometry(idx_pref).preferred_cl_res = 1;
                else % if all criteria fail, assume simple (as a catch)
                    idx_pref = strcmpi({obj.multi_allometry.type},'restrictive') & strcmpi({obj.multi_allometry.parameter},'CL_mLminkg') & strcmpi({obj.multi_allometry.method},'Simple');
                    obj.multi_allometry(idx_pref).preferred_cl_res = 1;
                end
                if round(nonres_exp,2) <= 0.70
                    idx_pref = strcmpi({obj.multi_allometry.type},'non-restrictive') & strcmpi({obj.multi_allometry.parameter},'CL_mLminkg') & strcmpi({obj.multi_allometry.method},'Simple');
                    obj.multi_allometry(idx_pref).preferred_cl_nonres = 1;
                elseif round(nonres_exp,2) > 0.70 && round(res_exp,2) <=1.00
                    idx_pref = strcmpi({obj.multi_allometry.type},'non-restrictive') & strcmpi({obj.multi_allometry.parameter},'CL_mLminkg') & strcmpi({obj.multi_allometry.method},'MLP');
                    obj.multi_allometry(idx_pref).preferred_cl_nonres = 1;
                elseif round(nonres_exp,2) > 1.00
                    idx_pref = strcmpi({obj.multi_allometry.type},'non-restrictive') & strcmpi({obj.multi_allometry.parameter},'CL_mLminkg') & strcmpi({obj.multi_allometry.method},'BrW');
                    obj.multi_allometry(idx_pref).preferred_cl_nonres = 1;
                else % if all criteria fail, assume simple (as a catch)
                    idx_pref = strcmpi({obj.multi_allometry.type},'non-restrictive') & strcmpi({obj.multi_allometry.parameter},'CL_mLminkg') & strcmpi({obj.multi_allometry.method},'Simple');
                    obj.multi_allometry(idx_pref).preferred_cl_nonres = 1;
                end
            end
            
        end
    end
end

% MLP (years)
% mouse = 2.7
% rat = 4.7
% guinea pig = 9
% marmoset = 10
% dog = 20
% monkey = 22
% minipig = 27
% human = 93
% rabbit = 8
