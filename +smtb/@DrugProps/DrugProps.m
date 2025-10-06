classdef DrugProps
    % DRUGPROPS A class that implements a drug properties object in the 
    % SMTB toolbox. A DrugProps object can contain basic 
    % drug properties such as:
    %    gcn                     %GSK compound number
    %    smiles                  %Smiles string
    %    mw                      %Molcular weight
    %    mf                      %Molecular formula
    %    pKa                     %pKas
    %    logP                    %logP
    %    bp                      %Blood:Plasma Ratio
    %    fup                     %Fraction unbound in plasma
    %    perm                    %Permeability
    %    met                     %In vitro metabolism
    %    pk_kate                 %PK Summary from KATE
    %
    % dprops = DrugProps(gcn), where gcn is a string containg 
    %   the GSK compound
    %   number or a cell array of multiple gsk compound numbers. This pulls
    %   the ADMET, ADAMANTIS, and KATE data typically relevent to PK
    %   modeling.
    %
    % dprops = DrugProps(smiles_string,'smiles'), automatically triggers a
    % search by smiles string rather than the DEX lookup service
    %
    % See also PARAMETER, INVIVOPROPS
    
    properties
        name            = '';
        gcn             = '';            %GSK compound number
        smiles          = '';            %Smiles string
        mw              = [];            %Molcular weight
        mf              = '';            %Molecular formula
        Deff            = struct('value',{}, ...
            'source',{});   % Effective diffusion coefficient
        sf              = []; % Solublization factor for ionized species
        pKa             = struct('acidic',{},'basic',{}, ...
            'cmpd_form',{},'source',{},'study',{},'lnb_ref',{}, ...
            'prep_lnb_ref',{},'protocol',{},'date_experiment',{}, ...
            'date_record',{},'date_updated',{},'comments',{}); % pKas
        logP            = struct('cmpd',{},'value',{},'source',{}, ...
            'lnb_ref',{},'prep_lnb_ref',{},'date_experiment',{}, ...
            'date_record',{},'date_updated',{}); % logP
        bp              = struct('cmpd',{},'value',{},'std',{},'n', ...
            {},'species',{},'details',{},'source',{},'flag',{}, ...
            'study',{},'lnb_ref',{},'prep_lnb_ref',{},'protocol',{}, ...
            'date_experiment',{},'date_record',{},'date_updated',{}, ...
            'comments',{}) % Blood:Plasma Ratio
        fup             = struct('cmpd',{},'value',{},'species',{}, ...
            'details',{},'source',{},'method',{},'study',{},'lnb_ref',{}, ...
            'prep_lnb_ref',{},'protocol',{},'date_experiment',{}, ...
            'date_record',{}, ...
            'date_updated',{},'comments',{}); % Fraction unbound in plasma
        bind            = struct('cmpd',{},'matrix',{},'fu',{}, ...
            'species',{},'details',{},'source',{},'method',{}, ...
            'study',{},'lnb_ref',{},'prep_lnb_ref',{},'protocol',{}, ...
            'date_experiment',{},'date_record',{},'date_updated',{}, ...
            'comments',{}); % Protein Binding
        perm            = struct('cmpd',{},'Modifier',{},'value',{}, ...
            'unit',{},'efflux_ratio',{},'pct_recovery',{},'type',{}, ...
            'category',{},'source',{},'study',{},'lnb_ref',{}, ...
            'prep_lnb_ref',{},'protocol',{},'date_experiment',{}, ...
            'date_record',{},'date_updated',{}, ...
            'comments',{}); % Permeability
        sol             = struct('cmpd',{},'matrix',{},'pH',{}, ...
            'Modifier',{},'value',{},'unit',{},'source',{},'lnb_ref', ...
            {},'prep_lnb_ref',{},'protocol',{},'date_experiment',{}, ...
            'date_record',{},'date_updated',{},'comments',{});% Solubility
        met             = struct('cmpd',{},'Modifier',{}, ...
            'CL',{},'unit',{},'species',{},'system',{},'details',{}, ...
            'source',{},'study',{},'lnb_ref',{},'prep_lnb_ref',{}, ...
            'protocol',{},'date_experiment',{},'date_record',{}, ...
            'date_updated',{},'comments',{});% In vitro metabolism
        enz             = struct('Enzyme',{},'Substrate', ...
            {},'Inhibitor',{},'source',{}); % CYP, UGT, and
        % PGP substrate/inhibition information
        descriptors     = struct('value',{},'type',{}, ...
            'source',{}); %molecular descriptors
        pk_kate         = struct('cmpd',{},'species',{}, ...
            'strain',{},'sex',{},'animal_id',{},'n',{},'dose_mgkg',{}, ...
            'route',{},'formulation',{},'matrix',{},'F_mod',{}, ...
            'F',{},'F_std',{},'Cmax_mod',{},'Cmax_ngml',{}, ...
            'Cmax_std',{},'Tmax_mod',{},'Tmax_hr',{},'Tmax_std',{}, ...
            'AUC_Inf_mod',{},'AUC_Inf_nghrmL',{},'AUC_Inf_std',{}, ...
            'AUC_t_mod',{},'AUC_t_nghrmL',{},'AUC_t_std',{}, ...
            'Clast_mod',{},'Clast_ngml',{},'Clast_std',{},'Tlast_mod',{}, ...
            'Tlast_hr',{},'TBdy_CL_mod',{},'TBdy_CL_mlmin_kg',{}, ...
            'TBdy_CL_std',{},'Vdss_mod',{},'Vdss_Lkg',{},'Vdss_std',{}, ...
            'T1_2_mod',{},'T1_2_hr',{},'T1_2_std',{},'MRT_mod',{}, ...
            'MRT_hr',{},'MRT_std',{},'FA_HPV_mod',{},'FA_HPV',{}, ...
            'FA_HPV_std',{},'project',{},'study',{},'study_design',{}, ...
            'lnb_ref',{},'prep_lnb_ref',{},'protocol_ID',{}, ...
            'pk_start_date',{},'pk_date_record',{},'pk_date_updated', ...
            {},'comments',{});% In vivo PK data from KATE
        pk_kate_summary = struct.empty;
        derived_props   = struct('type',{},'value',{},'units',{}, ...
            'source',{});
        fxn_status      = struct('service',{},'time_s',{}, ...
            'FromCache',{},'ReturnedValues',{},'comments',{});
    end
    
    methods
        function obj = DrugProps(name,varargin)
%             mtbc = MTBToolboxProxy;
            options = struct;
            options.searchType      = 'gcn';
            options.lookupGCN       = true;
            options.searchN         = 'single';
            options.getADMET        = true;
            options.getKATEinvitro  = true;
            options.getKATEpk       = true;
            options.getDPS          = true;
            options.addDPS          = false;
            options.addKATE         = false;
            options.addADMET        = false;
            options.PKoperation     = 'geometric mean';
            options.UseParallel     = false;
            
            for i = 1:2:length(varargin)
                switch strrep(strrep(lower(varargin{i}),'_',''),' ','')
                    case 'searchtype',      options.searchType      = smtb.useful.typecheck(varargin{i+1},'char','searchType');
                    case 'lookupgcn',       options.lookupGCN       = smtb.useful.typecheck(varargin{i+1},'logical','lookupGCN');
                    case 'searchn',         options.searchN         = smtb.useful.typecheck(varargin{i+1},'char','searchN');
                    case 'getadmet',        options.getADMET        = smtb.useful.typecheck(varargin{i+1},'logical','getADMET');
                    case 'getkateinvitro',  options.getKATEinvitro  = smtb.useful.typecheck(varargin{i+1},'logical','getKATEinvitro');
                    case 'getkatepk',       options.getKATEpk       = smtb.useful.typecheck(varargin{i+1},'logical','getKATEpk');
                    case 'getkate',         options.getKATEpk       = smtb.useful.typecheck(varargin{i+1},'logical','getKATE');
                                            options.getKATEinvitro  = smtb.useful.typecheck(varargin{i+1},'logical','getKATE');
                    case 'getdps',          options.getDPS          = smtb.useful.typecheck(varargin{i+1},'logical','getDPS');
                    case 'adddps',          options.addDPS          = smtb.useful.typecheck(varargin{i+1},'logical','addDPS');
                    case 'addkate',         options.addKATE         = smtb.useful.typecheck(varargin{i+1},'logical','addKATE');
                    case 'addadmet',        options.addADMET        = smtb.useful.typecheck(varargin{i+1},'logical','addADMET');
                    case 'pkoperation',     options.PKoperation     = smtb.useful.typecheck(varargin{i+1},'char','PK operation');
                    case 'useparallel',     warning('Option to use parallel computing for DrugProps is currently disabled.')
                    otherwise, error('Invalid input %s.  Valid inputs are getADMET, getKATE, getKATEinvitro, getKATEpk, getDPS, or SearchType.',varargin{i})
                end
            end
            
            if any([options.addDPS options.addKATE options.addADMET])
                if ~isa(name,'DrugProps')
                    error('Input must be a class DrugProps if seeking to add DPS, KATE, or ADMET properties to existing DrugProps.')
                else
                    obj = name;
                    if options.addDPS
                       try
                           obj = obj.getDPSProperties(obj); 
                       catch
                           warning('Unable to get GSK QSAR properties for %s',obj.gcn);
                           obj = obj.addPropertyValue('fxn_status','service','BuildMasterADMETCompoundProfile','ReturnedValues',0,'comments','Failure in DrugProps/getDPSProperties.');
                       end
                    else
                        mw = getSelectDPSmodelOutputs(strrep(obj.smiles,'#','%23'),'pis.mw');
                        mf = getSelectDPSmodelOutputs(strrep(obj.smiles,'#','%23'),'pis.mf');
                        if ~isempty(mw)
                            obj.mw = str2double(mw(strcmpi({mw.PropertyName},'pis.mw')).PropertyValue);
                        else
                            obj.mw = NaN;
                        end
                        if ~isempty(mf)
                            obj.mf = mf(strcmpi({mf.PropertyName},'pis.mf')).PropertyValue;
                        else
                            obj.mf = 'Error retrieving MF';
                        end
                    end
                    if options.addADMET
                       try
                            obj = obj.getADMETProperties(obj);
                       catch
                            warning('Unable to get ADMET properties.');
                            obj = obj.addPropertyValue('fxn_status','service','BuildMasterADMETCompoundProfile','ReturnedValues',0,'comments','Failure in DrugProps/getADMETProperties.');
                       end
                    end
                    if options.addKATE
                       try
                            obj = obj.getKATEPKProperties(obj);
                       catch
                           warning('Unable to get KATE PK results for %s',obj.gcn);
                       end
                       try
                           pk_bin = binKATEpkdata(obj,'operation',options.PKoperation);
                           obj.pk_kate_summary = rmfield(pk_bin,{'gcn','cmpd','upper_dose','lower_dose','pk_start_date','pk_date_record','pk_date_updated','comments','formulation'});
                       catch
                           warning('Error binning KATE PK data for %s',obj.gcn);
                       end
                       try
                            obj = obj.getKATEinVitroProperties(obj);
                       catch
                           warning('Unable to get KATE in vitro properties for %s',obj.gcn);
                       end
                    end
                end
                return
            end
            
            if iscell(name)
                obj = cell(1,length(name));
                for i = 1:length(name)
                    if logical(isGSKcn(name{i}))
                       try
                           obj{i} = DrugProps(name{i},'searchType','gcn','searchN','single',...
                               'getADMET',false,'getKATEpk',options.getKATEpk,'getKATEinvitro',options.getKATEinvitro,'getDPS',options.getDPS);
                       catch
                           warning('Unable to retrieve DrugProps for compound %s',name{i})
                       end
                    else
                        try
                            obj{i} = DrugProps(name{i},'searchType','smiles','searchN','single',...
                                'getADMET',false,'getKATEpk',options.getKATEpk,'getKATEinvitro',options.getKATEinvitro,'getDPS',options.getDPS);
                        catch
                            warning('Unable to retrieve DrugProps for compound %s',name{i})
                        end
                    end
                end
                obj = [obj{1:length(obj)}];
                if options.getADMET
                    try
                        obj = obj.getADMETProperties(obj); % Run ADMET service separate as a single batch. Switch statement prevents it from running single compounds in for loop above.
                    catch
                        warning('Unable to get ADMET properties using multi compound search');
                        for i = 1:length(obj)
                            obj(i) = obj(i).addPropertyValue('fxn_status','service','BuildMasterADMETCompoundProfile','ReturnedValues',0,'comments','Failed to return multi-compound ADMET Predictor profile.');
                        end
                    end
                end
                for i = 1:length(obj)
                    if ~any(contains({obj(i).perm.source},'QSAR Peff NN'))
                        try
                            obj(i) = obj(i).addPropertyValue('perm','cmpd',obj(i).name,'value',PeffNN(obj(i)),'unit','x10^-4 cm/s','type','Peff','source','QSAR Peff NN');
                        catch
                        end
                    end
                end
                return;
            end
           
           obj.name    = smtb.useful.typecheck(name,'char','name');
           
           switch lower(options.searchType)
               case 'gcn'
                   cl = smtb.webserv.searchCompound(name);
                   switch length(cl)
                       case 1
                           obj.gcn     = cl{:};
                           if ~strcmpi(obj.gcn,'not_found')
                               obj = obj.addPropertyValue('fxn_status','service','searchCompound','ReturnedValues',1);
                               try
                                   obj.smiles  = smtb.webserv.GetSMILESFromCLrest(cl{:});
                               catch
                                   obj.smiles = 'Error retrieving SMILES string';
                               end
                           else
                               obj = obj.addPropertyValue('fxn_status','service','searchCompound','ReturnedValues',0,'comments',sprintf('compound %s not found in DISH or CODS',obj.name));
                           end
                       case 0
                           fprintf('No compounds found in DEX or CODS for "%s". KATE data may not be available \n',name)
                           obj.gcn = name;
                           try
                               obj.smiles =  smtb.webserv.GetSMILESFromCLrest(name);
                           catch
                               obj.smiles = 'Error retrieving SMILES string';
                           end
                       otherwise
                           fprintf('Multiple compounds found in DEX lookup for "%s". Please choose from the following: \n',name)
                           for j = 1:length(cl), fprintf('\t%s\n',cl{j}); end
                           return
                   end
               case 'smiles'
                   obj.smiles = name;
                   obj.name = '';
                   if options.lookupGCN
                       gcn = getGCNfromSMILES(name);
                       obj.name = gcn;
                       obj.gcn = gcn;
                       if isempty(gcn)
                           warning('Compound number not found for SMILES string. KATE data will not be returned.');
                       end
                   end
               otherwise, error('Invalid search qualifier "%s"',searchType);
           end
                      
           switch lower(options.searchN)
               case 'single'
                   if ~isempty(obj.smiles)
                       if options.getADMET
                           try
                                obj = obj.getADMETProperties(obj);
                           catch
                                warning('Unable to get ADMET properties.');
                                obj = obj.addPropertyValue('fxn_status','service','BuildMasterADMETCompoundProfile','ReturnedValues',0,'comments','Failure in DrugProps/getADMETProperties.');
                           end
                       end
                   else
                       obj = obj.addPropertyValue('fxn_status','service','BuildMasterADMETCompoundProfile','ReturnedValues',0,'comments','empty SMILES string - ADMET Predictor properties not retrieved');
                   end
               otherwise
           end
           
           if options.getDPS
               if ~isempty(obj.smiles)
                   try
                       obj = obj.getDPSProperties(obj); 
                   catch
                       warning('Unable to get GSK QSAR properties for %s',obj.gcn);
                       obj = obj.addPropertyValue('fxn_status','service','BuildMasterADMETCompoundProfile','ReturnedValues',0,'comments','Failure in DrugProps/getDPSProperties.');
                   end
                   try
                       mw = smtb.webserv.getSelectDPSmodelOutputs(strrep(obj.smiles,'#','%23'),'pis.mw');
                       mf = smtb.webserv.getSelectDPSmodelOutputs(strrep(obj.smiles,'#','%23'),'pis.mf');
                       if ~isempty(mw)
                           obj.mw = str2double(mw(strcmpi({mw.PropertyName},'pis.mw')).PropertyValue);
                       else
                           obj.mw = NaN;
                       end
                       if ~isempty(mf)
                           obj.mf = mf(strcmpi({mf.PropertyName},'pis.mf')).PropertyValue;
                       else
                           obj.mf = 'Error retrieving MF';
                       end
                   catch
                       warning('Unable to retrieve MW and MF for %s',obj.gcn)
                   end
               else
                   obj = obj.addPropertyValue('fxn_status','service','BuildMasterADMETCompoundProfile','ReturnedValues',0,'comments','empty SMILES string - DPS properties not retrieved');
               end
           end
           
           if options.getKATEpk
               if ~isempty(obj.gcn) && ~strcmpi(obj.gcn,'not_found')
                   try
                        obj = obj.getKATEPKProperties(obj);
                   catch
                       warning('Unable to get KATE PK results for %s',obj.gcn);
                       obj = obj.addPropertyValue('fxn_status','service','GetInVivoPharmacokineticsResults','ReturnedValues',0,'comments','Failure in DrugProps/getKATEPKproperties.');
                   end
                   try
                       pk_bin = binKATEpkdata(obj,'operation',options.PKoperation);
                       obj.pk_kate_summary = rmfield(pk_bin,{'gcn','cmpd','upper_dose','lower_dose','pk_start_date','pk_date_record','pk_date_updated','comments','formulation'});
                   catch
                       warning('Error binning KATE PK data for %s',obj.gcn);
                   end
               else
                   obj = obj.addPropertyValue('fxn_status','service','GetInVivoPharmacokineticsResults','ReturnedValues',0,'comments',sprintf('compound %s not found in DISH or CODS',obj.name));
               end
           end
           
           if options.getKATEinvitro
               if ~isempty(obj.gcn) && ~strcmpi(obj.gcn,'not_found')
                   try
                        obj = obj.getKATEinVitroProperties(obj);
                   catch
                       warning('Unable to get KATE in vitro properties for %s',obj.gcn);
                       obj = obj.addPropertyValue('fxn_status','service','KATESearch','ReturnedValues',0,'comments','Failure in DrugProps/getKATEinVitroProperties.');
                   end
               else
                   obj = obj.addPropertyValue('fxn_status','service','KATESearch','ReturnedValues',0,'comments',sprintf('compound %s not found in DISH or CODS',obj.name));
               end
           end
           
           obj = obj.getDatamartProperties(obj);
           
           try
               obj = obj.addPropertyValue('perm','cmpd',obj.name,'value',PeffNN(obj),'unit','x10^-4 cm/s','type','Peff','source','QSAR Peff NN');
           catch
           end
               
           [~,idx] = sort({obj.fup.species}); obj.fup = obj.fup(idx);
           [~,idx] = sort({obj.bp.species});  obj.bp  = obj.bp(idx);
           [~,idx] = sort({obj.met.system});  obj.met = obj.met(idx);
           [~,idx] = sort({obj.met.species}); obj.met = obj.met(idx);
           
        end
        
        function obj = addPropertyValue(obj,prop,varargin)
            idx = length(obj.(prop))+1;
            for i = 1:2:length(varargin)
                obj.(prop)(idx).(varargin{i}) = varargin{i+1};
            end
        end
        
        function cl = getInVitroClearance(obj,varargin)
            cl = obj.met;
            for i = 1:2:length(varargin)
                switch lower(varargin{i})
                    case 'species', cl = cl(strcmpi({cl.species},varargin{i+1}));
                    case 'system',  cl = cl(strcmpi({cl.system},varargin{i+1}));
                    case 'source',  cl = cl(contains({cl.source},varargin{i+1},'IgnoreCase',true));
                    otherwise, error('Invalid input "%s"',varargin{i})
                end
            end
        end
        
        function [cl,cl_vals,matrix,strain,study_design,n,std,date_upd,dose] = getInVivoClearance(obj,species,varargin)
%             pk = obj.pk_kate_summary(strcmpi({obj.pk_kate_summary.species},species));
            pk = obj.pk_kate(strcmpi({obj.pk_kate.species},species));
            pk = pk([pk.dose_mgkg]>0);
            idx_iv = ismember(lower({pk.route}),{'iv' 'ivbolus' 'ivinf' 'ivinfusion' 'iv-inf'});
            pk = pk(idx_iv);
            matrix_opts = {'plasma','blood'};
            study_design_opts = {'discrete','cassette'};
            matrix = 'plasma';
            study_design = 'all';
            strain = '';
            operation = 'geometric mean';
            for i = 1:2:length(varargin)
                switch lower(varargin{i})
                    case 'matrix',      matrix = varargin{i+1};
                    case 'strain',      strain = varargin{i+1};
                    case 'study_design',study_design = varargin{i+1};
                    case 'operation',   operation = varargin{i+1};
                end
            end
            if ~isempty(pk)
                if strcmpi(study_design,'all')
                    idx_study_design = ones(1,length(pk));
                else
                    idx_study_design = strcmpi({pk.study_design},study_design);
                    if ~any(idx_study_design)
                        idx_study_design = strcmpi({pk.study_design},study_design_opts{~strcmpi(study_design_opts,study_design)});
                    end
                end
                idx_matrix = strcmpi({pk.matrix},matrix);
                if ~any(idx_matrix)
                    idx_matrix = strcmpi({pk.matrix},matrix_opts{~strcmpi(matrix_opts,matrix)});
                end
                if ~isempty(strain)
                    idx_strain = strcmpi({pk.strain},strain);
                    if ~any(idx_strain)
                        idx_strain = ones(1,length(pk));
                    end
                else
                    idx_strain = ones(1,length(pk));
                end
                idx_with_vals = cellfun(@(x) ~isempty(x),{pk.TBdy_CL_mlmin_kg},'UniformOutput',false);
                idx = idx_study_design & idx_matrix & idx_strain & logical(cell2mat(idx_with_vals));
                if any(idx)
%                     cl_vals = [pk(idx).CL_mLminkg];
                    cl_vals = [pk(idx).TBdy_CL_mlmin_kg];
                    matrix = {pk(idx).matrix};
                    strain = {pk(idx).strain};
                    study_design = {pk(idx).study_design};
                    n = [pk(idx).n];
                    std = [pk(idx).TBdy_CL_std];
                    date_upd = {pk(idx).pk_date_updated};
                    dose = [pk(idx).dose_mgkg];
                    switch strrep(strrep(strrep(lower(operation),'_',''),'-',''),' ','')
                        case 'geometricmean'
                            cl = geomean(cl_vals);
                        case {'arithmeticmean','average'}
                            cl = mean(cl_vals);
                        case 'median'
                            cl = median(cl_vals);
                        otherwise
                            error('%s is invalid operation option. Please choose from geometric mean, arithmetic mean, or median.')
                    end
                else
                    cl = [];
                    cl_vals = [];
                    matrix = '';
                    strain = '';
                    study_design = '';
                    n = [];
                    std = [];
                    date_upd = '';
                    dose = [];
                end
            else
                cl = [];
                cl_vals = [];
                matrix = '';
                strain = '';
                study_design = '';
                n = [];
                std = [];
                date_upd = '';
                dose = [];
            end
        end
        
        function [vss,vss_vals,matrix,strain,study_design,n,std,date_upd,dose] = getInVivoVss(obj,species,varargin)
%             pk = obj.pk_kate_summary(strcmpi({obj.pk_kate_summary.species},species));
            pk = obj.pk_kate(strcmpi({obj.pk_kate.species},species));
            pk = pk([pk.dose_mgkg]>0);
            idx_iv = ismember(lower({pk.route}),{'iv' 'ivbolus' 'ivinf' 'ivinfusion' 'iv-inf'});
            pk = pk(idx_iv);
            matrix_opts = {'plasma','blood'};
            study_design_opts = {'discrete','cassette'};
            matrix = 'plasma';
            study_design = 'all';
            strain = '';
            operation = 'geometric mean';
            for i = 1:2:length(varargin)
                switch lower(varargin{i})
                    case 'matrix',      matrix = varargin{i+1};
                    case 'strain',      strain = varargin{i+1};
                    case 'study_design',study_design = varargin{i+1};
                    case 'operation',   operation = varargin{i+1};
                end
            end
            if ~isempty(pk)
                if strcmpi(study_design,'all')
                    idx_study_design = ones(1,length(pk));
                else
                    idx_study_design = strcmpi({pk.study_design},study_design);
                    if ~any(idx_study_design)
                        idx_study_design = strcmpi({pk.study_design},study_design_opts{~strcmpi(study_design_opts,study_design)});
                    end
                end
                idx_matrix = strcmpi({pk.matrix},matrix);
                if ~any(idx_matrix)
                    idx_matrix = strcmpi({pk.matrix},matrix_opts{~strcmpi(matrix_opts,matrix)});
                end
                if ~isempty(strain)
                    idx_strain = strcmpi({pk.strain},strain);
                    if ~any(idx_strain)
                        idx_strain = ones(1,length(pk));
                    end
                else
                    idx_strain = ones(1,length(pk));
                end
                idx_with_vals = cellfun(@(x) ~isempty(x),{pk.Vdss_Lkg},'UniformOutput',false);
                idx = idx_study_design & idx_matrix & idx_strain & logical(cell2mat(idx_with_vals));
                if any(idx)
%                     cl_vals = [pk(idx).CL_mLminkg];
                    vss_vals = [pk(idx).Vdss_Lkg];
                    matrix = {pk(idx).matrix};
                    strain = {pk(idx).strain};
                    study_design = {pk(idx).study_design};
                    n = [pk(idx).n];
                    std = [pk(idx).Vdss_std];
                    date_upd = {pk(idx).pk_date_updated};
                    dose = [pk(idx).dose_mgkg];
                    switch strrep(strrep(strrep(lower(operation),'_',''),'-',''),' ','')
                        case 'geometricmean'
                            vss = geomean(vss_vals);
                        case {'arithmeticmean','average'}
                            vss = mean(vss_vals);
                        case 'median'
                            vss = median(vss_vals);
                        otherwise
                            error('%s is invalid operation option. Please choose from geometric mean, arithmetic mean, or median.')
                    end
                else
                    vss = [];
                    vss_vals = [];
                    matrix = '';
                    strain = '';
                    study_design = '';
                    n = [];
                    std = [];
                    date_upd = '';
                    dose = [];
                end
            else
                vss = [];
                vss_vals = [];
                matrix = '';
                strain = '';
                study_design = '';
                n = [];
                std = [];
                date_upd = '';
                dose = [];
            end
        end
        
        function [Thalf,Thalf_vals,matrix,strain,study_design,n,std,date_upd,dose] = getInVivoThalf(obj,species,varargin)
%             pk = obj.pk_kate_summary(strcmpi({obj.pk_kate_summary.species},species));
            pk = obj.pk_kate(strcmpi({obj.pk_kate.species},species));
            pk = pk([pk.dose_mgkg]>0);
            idx_iv = ismember(lower({pk.route}),{'iv' 'ivbolus' 'ivinf' 'ivinfusion' 'iv-inf'});
            pk = pk(idx_iv);
            matrix_opts = {'plasma','blood'};
            study_design_opts = {'discrete','cassette'};
            matrix = 'plasma';
            study_design = 'all';
            strain = '';
            operation = 'geometric mean';
            for i = 1:2:length(varargin)
                switch lower(varargin{i})
                    case 'matrix',      matrix = varargin{i+1};
                    case 'strain',      strain = varargin{i+1};
                    case 'study_design',study_design = varargin{i+1};
                    case 'operation',   operation = varargin{i+1};
                end
            end
            if ~isempty(pk)
                if strcmpi(study_design,'all')
                    idx_study_design = ones(1,length(pk));
                else
                    idx_study_design = strcmpi({pk.study_design},study_design);
                    if ~any(idx_study_design)
                        idx_study_design = strcmpi({pk.study_design},study_design_opts{~strcmpi(study_design_opts,study_design)});
                    end
                end
                idx_matrix = strcmpi({pk.matrix},matrix);
                if ~any(idx_matrix)
                    idx_matrix = strcmpi({pk.matrix},matrix_opts{~strcmpi(matrix_opts,matrix)});
                end
                if ~isempty(strain)
                    idx_strain = strcmpi({pk.strain},strain);
                    if ~any(idx_strain)
                        idx_strain = ones(1,length(pk));
                    end
                else
                    idx_strain = ones(1,length(pk));
                end
                idx_with_vals = cellfun(@(x) ~isempty(x),{pk.T1_2_hr},'UniformOutput',false);
                idx = idx_study_design & idx_matrix & idx_strain & logical(cell2mat(idx_with_vals));
                if any(idx)
%                     cl_vals = [pk(idx).CL_mLminkg];
                    Thalf_vals = [pk(idx).T1_2_hr];
                    matrix = {pk(idx).matrix};
                    strain = {pk(idx).strain};
                    study_design = {pk(idx).study_design};
                    n = [pk(idx).n];
                    std = [pk(idx).T1_2_std];
                    date_upd = {pk(idx).pk_date_updated};
                    dose = [pk(idx).dose_mgkg];
                    switch strrep(strrep(strrep(lower(operation),'_',''),'-',''),' ','')
                        case 'geometricmean'
                            Thalf = geomean(Thalf_vals);
                        case {'arithmeticmean','average'}
                            Thalf = mean(Thalf_vals);
                        case 'median'
                            Thalf = median(Thalf_vals);
                        otherwise
                            error('%s is invalid operation option. Please choose from geometric mean, arithmetic mean, or median.')
                    end
                else
                    Thalf = [];
                    Thalf_vals = [];
                    matrix = '';
                    strain = '';
                    study_design = '';
                    n = [];
                    std = [];
                    date_upd = '';
                    dose = [];
                end
            else
                Thalf = [];
                Thalf_vals = [];
                matrix = '';
                strain = '';
                study_design = '';
                n = [];
                std = [];
                date_upd = '';
                dose = [];
            end
        end
        
        function [mrt,mrt_vals,matrix,strain,study_design,n,std,date_upd,dose] = getInVivoMRT(obj,species,varargin)
%             pk = obj.pk_kate_summary(strcmpi({obj.pk_kate_summary.species},species));
            pk = obj.pk_kate(strcmpi({obj.pk_kate.species},species));
            pk = pk([pk.dose_mgkg]>0);
            idx_iv = ismember(lower({pk.route}),{'iv' 'ivbolus' 'ivinf' 'ivinfusion' 'iv-inf'});
            pk = pk(idx_iv);
            matrix_opts = {'plasma','blood'};
            study_design_opts = {'discrete','cassette'};
            matrix = 'plasma';
            study_design = 'all';
            strain = '';
            operation = 'geometric mean';
            for i = 1:2:length(varargin)
                switch lower(varargin{i})
                    case 'matrix',      matrix = varargin{i+1};
                    case 'strain',      strain = varargin{i+1};
                    case 'study_design',study_design = varargin{i+1};
                    case 'operation',   operation = varargin{i+1};
                end
            end
            if ~isempty(pk)
                if strcmpi(study_design,'all')
                    idx_study_design = ones(1,length(pk));
                else
                    idx_study_design = strcmpi({pk.study_design},study_design);
                    if ~any(idx_study_design)
                        idx_study_design = strcmpi({pk.study_design},study_design_opts{~strcmpi(study_design_opts,study_design)});
                    end
                end
                idx_matrix = strcmpi({pk.matrix},matrix);
                if ~any(idx_matrix)
                    idx_matrix = strcmpi({pk.matrix},matrix_opts{~strcmpi(matrix_opts,matrix)});
                end
                if ~isempty(strain)
                    idx_strain = strcmpi({pk.strain},strain);
                    if ~any(idx_strain)
                        idx_strain = ones(1,length(pk));
                    end
                else
                    idx_strain = ones(1,length(pk));
                end
                idx_with_vals = cellfun(@(x) ~isempty(x),{pk.MRT_hr},'UniformOutput',false);
                idx = idx_study_design & idx_matrix & idx_strain & logical(cell2mat(idx_with_vals));
                if any(idx)
%                     cl_vals = [pk(idx).CL_mLminkg];
                    mrt_vals = [pk(idx).MRT_hr];
                    matrix = {pk(idx).matrix};
                    strain = {pk(idx).strain};
                    study_design = {pk(idx).study_design};
                    n = [pk(idx).n];
                    std = [pk(idx).MRT_std];
                    date_upd = {pk(idx).pk_date_updated};
                    dose = [pk(idx).dose_mgkg];
                    switch strrep(strrep(strrep(lower(operation),'_',''),'-',''),' ','')
                        case 'geometricmean'
                            mrt = geomean(mrt_vals);
                        case {'arithmeticmean','average'}
                            mrt = mean(mrt_vals);
                        case 'median'
                            mrt = median(mrt_vals);
                        otherwise
                            error('%s is invalid operation option. Please choose from geometric mean, arithmetic mean, or median.')
                    end
                else
                    mrt = [];
                    mrt_vals = [];
                    matrix = '';
                    strain = '';
                    study_design = '';
                    n = [];
                    std = [];
                    date_upd = '';
                    dose = [];
                end
            else
                mrt = [];
                mrt_vals = [];
                matrix = '';
                strain = '';
                study_design = '';
                n = [];
                std = [];
                date_upd = '';
                dose = [];
            end
        end
        
        function [F,F_vals,matrix,strain,study_design,n,std,date_upd,dose] = getInVivoF(obj,species,varargin)
%             pk = obj.pk_kate_summary(strcmpi({obj.pk_kate_summary.species},species));
            pk = obj.pk_kate(strcmpi({obj.pk_kate.species},species));
            pk = pk([pk.dose_mgkg]>0);
            idx_iv = ismember(lower({pk.route}),{'po','oral','po(capsule)','po(fast)','po(fed)','po(solution)','po(suspension)'});
            pk = pk(idx_iv);
            matrix_opts = {'plasma','blood'};
            study_design_opts = {'discrete','cassette'};
            matrix = 'plasma';
            study_design = 'all';
            strain = '';
            operation = 'geometric mean';
            for i = 1:2:length(varargin)
                switch lower(varargin{i})
                    case 'matrix',      matrix = varargin{i+1};
                    case 'strain',      strain = varargin{i+1};
                    case 'study_design',study_design = varargin{i+1};
                    case 'operation',   operation = varargin{i+1};
                end
            end
            if ~isempty(pk)
                if strcmpi(study_design,'all')
                    idx_study_design = ones(1,length(pk));
                else
                    idx_study_design = strcmpi({pk.study_design},study_design);
                    if ~any(idx_study_design)
                        idx_study_design = strcmpi({pk.study_design},study_design_opts{~strcmpi(study_design_opts,study_design)});
                    end
                end
                idx_matrix = strcmpi({pk.matrix},matrix);
                if ~any(idx_matrix)
                    idx_matrix = strcmpi({pk.matrix},matrix_opts{~strcmpi(matrix_opts,matrix)});
                end
                if ~isempty(strain)
                    idx_strain = strcmpi({pk.strain},strain);
                    if ~any(idx_strain)
                        idx_strain = ones(1,length(pk));
                    end
                else
                    idx_strain = ones(1,length(pk));
                end
                idx_with_vals = cellfun(@(x) ~isempty(x),{pk.F},'UniformOutput',false);
                idx = idx_study_design & idx_matrix & idx_strain & logical(cell2mat(idx_with_vals));
                if any(idx)
%                     cl_vals = [pk(idx).CL_mLminkg];
                    F_vals = [pk(idx).F];
                    matrix = {pk(idx).matrix};
                    strain = {pk(idx).strain};
                    study_design = {pk(idx).study_design};
                    n = [pk(idx).n];
                    std = [pk(idx).F_std];
                    date_upd = {pk(idx).pk_date_updated};
                    dose = [pk(idx).dose_mgkg];
                    switch strrep(strrep(strrep(lower(operation),'_',''),'-',''),' ','')
                        case 'geometricmean'
                            F = geomean(F_vals);
                        case {'arithmeticmean','average'}
                            F = mean(F_vals);
                        case 'median'
                            F = median(F_vals);
                        otherwise
                            error('%s is invalid operation option. Please choose from geometric mean, arithmetic mean, or median.')
                    end
                else
                    F = [];
                    F_vals = [];
                    matrix = '';
                    strain = '';
                    study_design = '';
                    n = [];
                    std = [];
                    date_upd = '';
                    dose = [];
                end
            else
                F = [];
                F_vals = [];
                matrix = '';
                strain = '';
                study_design = '';
                n = [];
                std = [];
                date_upd = '';
                dose = [];
            end
        end
        
        function pk = getPKdata(obj,varargin)
            pk = obj.pk_kate;
            sp = 'all';
            matrix = 'all';
            study_design = 'all';
            route = 'all';
            for i = 1:2:length(varargin)
                switch lower(varargin{i})
                    case 'matrix',      matrix = varargin{i+1};
                    case 'species',     sp = varargin{i+1};
                    case 'study_design',study_design = varargin{i+1};
                    case 'route',       route = varargin{i+1};
                end
            end
            
            if strcmpi(sp,'all')
                idx_sp = true(1,length(pk));
            else
                idx_sp = strcmpi({pk.species},sp);
            end
            
            if strcmpi(matrix,'all')
                idx_mat = true(1,length(pk));
            else
                idx_mat = strcmpi({pk.matrix},matrix);
            end
            
            if strcmpi(study_design,'all')
                idx_des = true(1,length(pk));
            else
                idx_des = strcmpi({pk.study_design},study_design);
            end
            
            switch lower(strrep(strrep(strrep(route,'-',''),'_',''),' ',''))
                case 'all'
                    idx_rt = true(1,length(pk));
                case {'po','oral'}
                    idx_rt = ismember(lower({pk.route}),{'po','oral','po(capsule)','po(fast)','po(fed)','po(solution)','po(suspension)'});
                case {'iv','ivinf','intravenous'}
                    idx_rt = ismember(lower({pk.route}),{'iv' 'ivbolus' 'ivinf' 'ivinfusion' 'iv-inf'});
                otherwise
                    idx_rt = true(length(pk),1);
            end
            
            idx = idx_sp & idx_mat & idx_des & idx_rt;
            
            pk = pk(idx);
            
        end
        
        function bp = getBP(obj,varargin)
            bp = obj.bp;
            for i = 1:2:length(varargin)
                switch lower(varargin{i})
                    case 'species', bp = bp(strcmpi({bp.species},varargin{i+1}));
                    case 'source',  bp = bp(contains({bp.source},varargin{i+1},'IgnoreCase',true));
                end
            end
        end
        
        function Fup = getFup(obj,varargin)
            Fup = obj.fup;
            for i = 1:2:length(varargin)
                switch lower(varargin{i})
                    case 'species', Fup = Fup(strcmpi({Fup.species}, ...
                            varargin{i+1}));
                    case 'source',  Fup = Fup(contains({Fup.source}, ...
                            varargin{i+1},'IgnoreCase',true));
                end
            end
        end
        
        function Fub = getFub(obj,varargin)
            Fub = obj.bind(strcmpi({obj.bind.matrix},'blood'));
            if ~isempty(Fub)
                for i = 1:2:length(varargin)
                   switch lower(varargin{i})
                       case 'species',  Fub = Fub(strcmpi({Fub.species},varargin{i+1}));
                       case 'source',   Fub = Fub(contains({Fub.source},varargin{i+1},'IgnoreCase',true));
                   end
                end
            end
        end
        
        function [sol,pH,source] = getSol(obj,varargin)
            sol_in = obj.sol;
            if ~isempty(sol_in)
                for i = 1:2:length(varargin)
                    switch lower(varargin{i})
                        case 'source',  sol_in = sol_in(contains({sol_in.source},varargin{i+1},'IgnoreCase',true));
                        case 'matrix',  sol_in = sol_in(strcmpi({sol_in.matrix},varargin{i+1}));
                    end
                end
            end
            sol = [sol_in.value];
            pH = [sol_in.pH];
            if isempty(pH), pH = -1 * ones(length(sol),1); end
            source = unique({sol_in.source});
            source = source{1};
        end
        
        function sol = calcpHSol(obj,pH)
            if length(pH) > 1
                sol = zeros(length(pH),1);
                for i = 1:length(pH), sol(i) = calcpHSol(obj,pH(i)); end
                return; 
            end
            iSol = obj.sol(strcmp({obj.sol.matrix},'Intrinsic')).value;
            [~,details]   = smtb.pbpkionfrac(getPreferredValue(obj,'pKa_a'),getPreferredValue(obj,'pKa_b'),pH);
            sol  = [iSol/details.F_all(1) iSol * obj.sf.value./details.F_all(2:end)];
            sol  = min(sol);
        end
        
        function disp(obj)
            for i = 1:length(obj)
                displayBasicInfo(obj(i));
                displayPhysChem(obj(i));
                displayBloodProps(obj(i));
                displayMetabolism(obj(i));
                displayClearanceRisk(obj(i));
                fprintf('PK Data Found in KATE, Use displayPKKate() for summary \n')
                if length(obj) > 1, fprintf('\n\n'); end
            end
        end
        
        function displayBasicInfo(obj)
            fprintf('Name:    %s \n',obj.name)
            fprintf('GSK #:   %s \n',obj.gcn)
            fprintf('SMILES:  %s \n',obj.smiles)
            fprintf('Formula: %s \n',obj.mf)
            fprintf('Mol Wgt: %0.2f Da \n',obj.mw)
        end
        
        function displayPhysChem(obj)
            if isempty (obj.pKa)
                fprintf('pKa:        (no data) \n')
            else
                fprintf('Acidic pKa: ');
                for i = 1:length(obj.pKa)
                    pKa_a = sprintf('%5.2f, ',obj.pKa(i).acidic); pKa_a = pKa_a(1:length(pKa_a)-2); 
                    if isempty(pKa_a), pKa_a = '(none)'; end
                    if i ==1, fprintf('%s (%s) \n',pKa_a,obj.pKa(i).source);
                    else,     fprintf('\t\t    %s (%s) \n',pKa_a,obj.pKa(i).source);
                    end
                end

                fprintf('Basic pKa: ');
                for i = 1:length(obj.pKa)
                    pKa_b = sprintf('%0.2f, ',obj.pKa(i).basic); pKa_b = pKa_b(1:length(pKa_b)-2); 
                    if isempty(pKa_b), pKa_b = '(none)'; end
                    if i ==1, fprintf('%s (%s) \n',pKa_b,obj.pKa(i).source);
                    else,     fprintf('\t\t   %s (%s) \n',pKa_b,obj.pKa(i).source);
                    end
                end
            end
            
            if isempty(obj.logP)
                fprintf('LogP:       (no data) \n');
            else
                fprintf('LogP:       %0.2f (%s)\n',obj.logP(1).value,obj.logP(1).source);
                for i = 2:length(obj.logP), fprintf('            %0.2f (%s)\n',obj.logP(i).value,obj.logP(i).source); end
            end

            fprintf('Solubility')
            if length(obj.sol) >= 1
                fprintf('\n\t Matrix \t  pH \t    Value \t\t\t Source \n')
                for i = 1:length(obj.sol)
                    if isempty(obj.sol(i).pH), fprintf('\t %-8s \t \t\t %6.3f %-8s \t %s\n',obj.sol(i).matrix,obj.sol(i).value,obj.sol(i).unit,obj.sol(i).source)
                    else,                      fprintf('\t %-8s \t %4.2f \t %6.3f %-8s \t %s\n',obj.sol(i).matrix,obj.sol(i).pH,obj.sol(i).value,obj.sol(i).unit,obj.sol(i).source)
                    end
                end
            else
                fprintf(':      (no data) \n');
            end
            
            fprintf('Permeability \n')
            for i = 1:length(obj.perm)
                fprintf('     %-45.45s: %2s%6.2f %s (%s,%s) \n',obj.perm(i).type,obj.perm(i).Modifier,obj.perm(i).value,obj.perm(i).unit,obj.perm(i).source,obj.perm(i).study);
            end
            
            fprintf('Efflux Ratio \n')
            idx = find(~contains({obj.perm.type},'Inh')&~contains({obj.perm.type},'BL/AP'));
            for i = 1:length(idx)
                if ~isempty(obj.perm(idx(i)).efflux_ratio)
                    fprintf('     %-35.35s: %6.2f (%s,%s) \n',obj.perm(idx(i)).type,obj.perm(idx(i)).efflux_ratio,obj.perm(i).source,obj.perm(i).study);
                end
            end
        end
        
        function displayBloodProps(obj)
            fprintf('Fraction unbound in plasma (fup) \n')
            for i = 1:length(obj.fup)
                fprintf('     %-6.6s: %5.2f%% (%s)\n',obj.fup(i).species,obj.fup(i).value.*100,obj.fup(i).source);
            end
            
            fprintf('Blood:Plasma Ratio')
            if ~isempty(obj.bp)
                fprintf('\n')
                for i = 1:length(obj.bp)
                    fprintf('     %-8.8s: %2.2f (%s) %s \n',obj.bp(i).species,obj.bp(i).value,obj.bp(i).source,obj.bp(i).flag);
                end
            else
                fprintf(':       (no data) \n');
            end
            
        end
        
        function displayMetabolism(obj)
            fprintf('In Vitro Metabolism')
            if ~isempty(obj.met)
                fprintf('\n')
                sp  = unique({obj.met.species});
                sys = unique({obj.met.system});
                src = unique({obj.met.source});
                for i = 1:length(sp)
                    for j = 1:length(sys)
                        for k = 1:length(src)
                            CL = [obj.met(strcmp({obj.met.species},sp{i}) & strcmp({obj.met.system},sys{j}) & strcmp({obj.met.source},src{k})).CL];
                            un = unique({obj.met(strcmp({obj.met.species},sp{i}) & strcmp({obj.met.system},sys{j})).unit});
%                             sc = unique({obj.met(strcmp({obj.met.species},sp{i}) & strcmp({obj.met.system},sys{j})).source});
                            sc = src{k};
                            switch length(CL)
                                case 0  % Do nothing, no data
                                case 1,    fprintf('     %-20.20s: %2.2f %s (%s) \n',[sp{i} ' ' sys{j}],10.^mean(log10(CL)),un{:},sc);
                                otherwise, fprintf('     %-20.20s: %2.2f %s (%2.2f-%2.2f %s) (%s) \n',[sp{i} ' ' sys{j}],10.^mean(log10(CL)),un{:},min(CL),max(CL),un{:},sc);
                            end
                        end
                    end
                end 
            else
                fprintf(':      (no data) \n');
            end
        end
        
        function displayClearanceRisk(obj)
            fprintf('Clearance Classification \n')
            if ~isempty(obj.derived_props)
                for i = 1:length(obj.derived_props)
                   if ~isnumeric(obj.derived_props(i).value)
                       fprintf('     %-30.30s: %s (%s) \n',obj.derived_props(i).type,obj.derived_props(i).value,obj.derived_props(i).source);
                   end
                end
            end
        end
        
        function displayPKKate(obj)
            fprintf('In Vivo PK')
            if ~isempty(obj.pk_kate)
                fprintf('\n\t');
                out = cell(2,13);
                out(1,:) = {'Compound   ',' Species ',' Dose ',' Route ',' Matrix ','  %F(%CV)  ','   Cmax(%CV) ','   Tmax(%CV) ','   AUC_inf(%CV) ','   Vd_ss(%CV) ','   CL(%CV) ','  Study    ','  Comments     '};
                out(2,:) = {'','',' mg/kg','','','',' ng/mL',' hr',' ng*hr/mL',' L/kg',' mL/min/kg','',''};
                
                for i = 1:length(obj.pk_kate)
                    d = obj.pk_kate(i);
                    out(2+i,:) = {d.cmpd,d.species,sprintf('%4.1f',d.dose_mgkg),d.route,d.matrix,'','','','','','','',''};
                    if ~isempty(d.F),               out{2+i,6}  = sprintf('%0.1f(%0.0f)',d.F,d.F_std/d.F*100); end
                    if ~isempty(d.Cmax_ngml),       out{2+i,7}  = sprintf('%0.0f(%0.0f)',d.Cmax_ngml,d.Cmax_sd/d.Cmax_ngml*100); end
                    if ~isempty(d.Tmax_hr),         out{2+i,8}  = sprintf('%0.1f(%0.0f)',d.Tmax_hr,d.Tmax_sd/d.Tmax_hr*100);end
                    if ~isempty(d.AUC_Inf_nghrmL),  out{2+i,9}  = sprintf('%0.0f(%0.0f)',d.AUC_Inf_nghrmL,d.AUC_Inf_std/d.AUC_Inf_nghrmL*100);end
                    if ~isempty(d.Vdss_Lkg),        out{2+i,10} = sprintf('%0.2f(%0.0f)',d.Vdss_Lkg,d.Vdss_std/d.Vdss_Lkg*100);end
                    if ~isempty(d.TBdy_CL_mlmin_kg),out{2+i,11} = sprintf('%0.1f(%0.0f)',d.TBdy_CL_mlmin_kg,d.TBdy_CL_std/d.TBdy_CL_mlmin_kg*100);end
                    if ~isempty(d.study),out{2+i,12} = sprintf(d.study);end
                    if ~isempty(d.comments),out{2+i,13} = sprintf(d.comments);end
                end
                for i = 1:size(out,1)
                    for k = 1:size(out,2)
                        if ~isempty(out{i,k}), output_parts = textscan(out{i,k},'%f(%f)'); else, output_parts = cell(2,1); end
                        nchar = max(cellfun(@length,out(:,k)));
                        if i < 3 || k < 6
                            tr_sp = blanks(floor((nchar - length(out{i,k}))/2)+1);
                            st_sp = blanks(ceil((nchar - length(out{i,k}))/2)+2);
                        else
                            cv_length = length(num2str(output_parts{2}));
                            tr_sp = blanks(3 - cv_length );
                            st_sp = blanks(nchar + 3 - length(out{i,k}) - length(tr_sp));
                        end
                        if (length(out{i,k}) + length(tr_sp) + length(st_sp)) ~= (nchar + 3), keyboard; end
                        fprintf('%s%s%s',st_sp,out{i,k},tr_sp);
                    end
                    fprintf('\n\t');
                end
            else
                fprintf('\n(no data) \n');
            end
        end
        
        %Defined elswhere
        res = allometricScaling(obj,weighting)
        res = summarizeF(obj);

    end
    
    methods (Static,Hidden,Access=private)
        function value = getPISproperties(obj,mtbc,prms)
            switch(prms)
                case 'mw'
                    try
                        prms = 'pis.mw';
                        dpsCalc = GetSelectCalculatedDerivedPropertiesDPS(mtbc, prms, obj.smiles);
                        value = str2double(dpsCalc.CmpdProperties(2).PropertyValue);
                    catch
                        value = NaN;
                        warning('Unable to retrieve molecular weight for "%s" from property information services',obj.gcn)
                    end
                case 'mf'
                    try 
                        prms = 'pis.mf';
                        dpsCalc = GetSelectCalculatedDerivedPropertiesDPS(mtbc, prms, obj.smiles);
                        value = dpsCalc.CmpdProperties(2).PropertyValue;
                    catch
                        value = '';
                        warning('Unable to retrieve molecular weight for "%s" from property information services',obj.gcn)
                    end
                      
                otherwise, error('Unrecognized PIS property')
            end
        end %getPISproperties
        
        function obj = getADMETProperties(obj)
            smiles = {obj(~contains({obj.smiles},'error','IgnoreCase',true)&~strcmpi({obj.smiles},'*')).smiles}'; %Remove cases where there was an error with SMILES or a * was returned as SMILES
            gcn = {obj(~contains({obj.smiles},'error','IgnoreCase',true)&~strcmpi({obj.smiles},'*')).gcn}';
            
            timer = tic;
            for i = 1:ceil(length(smiles)/100)
                if i==ceil(length(smiles)/100)
                    try
                        admet(((i*100)-99):length(smiles)) = smtb.webserv.getADMETProps(smiles(((i*100)-99):length(smiles)),gcn(((i*100)-99):length(smiles)));
                    catch
                        num = [(i*100)-99,length(smiles)];
                        warning('Unable to get ADMET properties using multi compound search for compounds %d through %d',num(1),num(2));
                    end
                else
                    try
                        admet(((i*100)-99):i*100) = getADMETProps(smiles(((i*100)-99):i*100),gcn(((i*100)-99):i*100));
                    catch
                        num = [(i*100)-99,length(obj)];
                        warning('Unable to get ADMET properties using multi compound search for compounds %d through %d',num(1),num(2));
                    end
                end
            end
            timer = toc(timer);
            
            for i = 1:length(admet)
                if ~isempty(admet(i).CmpdGeneral)
                    p = admet(i).CmpdPhysChem;
                    g = admet(i).CmpdGeneral;
                    m = admet(i).CmpdMetabolism;
                    t = admet(i).CmpdTopological;
                    c = admet(i).FromCache;
                    source = sprintf('QSAR ADMET v%2.1f',g.ADMETVersion);
                    
                    idx = find(strcmp({obj.smiles},g.Structure));
                    if isempty(idx),idx = find(strcmpi({obj.gcn},g.CmpdNumber)); end
                    
                    % Phys chem properties
                    if ~isempty(idx)
                        for j = 1:length(idx) % In cases where more than 1 compound is requested an
                            if isempty(obj(idx(j)).sf)
                                if ~any(strcmpi({obj(idx(j)).smiles},g.Structure)) && ~isempty(g.Structure)  % Replace salt SMILES with parent SMILES
                                    obj(idx(j)).smiles = g.Structure;
                                end
                                if ~isempty(g.Structure)
                                    obj(idx(j)) = obj(idx(j)).addPropertyValue('fxn_status','service','BuildMasterADMETCompoundProfile','time_s',timer,'FromCache',c,'ReturnedValues',1,'comments',admet(i).Comments);
                                    obj(idx(j)) = obj(idx(j)).addPropertyValue('logP','cmpd',g.CmpdNumber,'value',p.MlogP,'source',sprintf('%s Moriguchi',source));
                                    obj(idx(j)) = obj(idx(j)).addPropertyValue('logP','cmpd',g.CmpdNumber,'value',p.SlogP,'source',sprintf('%s S+logP',source));
                                    pKa_a = sscanf(p.SAcidic_pKa,'%f;'); pKa_a = sort(pKa_a);
                                    pKa_b = sscanf(p.SBasic_pKa,'%f;');  pKa_b = sort(pKa_b);
                                    obj(idx(j)) = obj(idx(j)).addPropertyValue('pKa','acidic',pKa_a,'basic',pKa_b,'source',source);
                                    obj(idx(j)) = obj(idx(j)).addPropertyValue('Deff','value', p.DiffCoef./1E5,'source',source);
                                    obj(idx(j)) = obj(idx(j)).addPropertyValue('descriptors','value',t.T_HydroR,'source',source,'type','r_he');
                                    obj(idx(j)) = obj(idx(j)).addPropertyValue('descriptors','value',t.T_Rads,'source',source,'type','r_s');
                                    
                                    eccs_class = '';
                                    switch strrep(strrep(lower(p.ECCS_Class),'_',''),' ','')
                                        case {'class1a','class2'},  eccs_class = 'Metabolism';
                                        case 'class1b',             eccs_class = 'Hepatic_Uptake';
                                        case {'class3a','class4'},  eccs_class = 'Renal';
                                        case 'class3b',             eccs_class = 'Heaptic_Uptake_or_Renal';
                                    end
                                    obj(idx(j)) = obj(idx(j)).addPropertyValue('derived_props','type','ECCS clearance','value',sprintf('%s_%s',p.ECCS_Class,eccs_class),'units','NA','source',sprintf('%s ECCS',source));

                                    % Solubility
                                    obj(idx(j)) = obj(idx(j)).addPropertyValue('sol','cmpd',g.CmpdNumber,'matrix','FaSSGF',   'pH',1.6,   'value',p.SFaSSGF,'unit','mg/mL','source',source);
                                    obj(idx(j)) = obj(idx(j)).addPropertyValue('sol','cmpd',g.CmpdNumber,'matrix','FaSSIF',   'pH',6.5,   'value',p.SFaSSIF,'unit','mg/mL','source',source);
                                    obj(idx(j)) = obj(idx(j)).addPropertyValue('sol','cmpd',g.CmpdNumber,'matrix','FeSSIF',   'pH',5.0,   'value',p.SFeSSIF,'unit','mg/mL','source',source);
                                    obj(idx(j)) = obj(idx(j)).addPropertyValue('sol','cmpd',g.CmpdNumber,'matrix','Aqueous',  'pH',p.SpH, 'value',p.SSw,    'unit','mg/mL','source',source);
                                    obj(idx(j)) = obj(idx(j)).addPropertyValue('sol','cmpd',g.CmpdNumber,'matrix','Intrinsic','pH',[],    'value',p.SIS,    'unit','mg/mL','source',source);

                                    % Ion solubilization factor
                                    obj(idx(j)).sf.value  = p.SSF;
                                    obj(idx(j)).sf.source = source;

                                    % Permeability
                                    obj(idx(j)) = obj(idx(j)).addPropertyValue('perm','cmpd',g.CmpdNumber,'value',p.SPeff,'unit','x10^-4 cm/s','type','Peff','category','Human Jejunum','source',source);
                                    obj(idx(j)) = obj(idx(j)).addPropertyValue('perm','cmpd',g.CmpdNumber,'value',p.SMDCK/10, 'unit','nm/s','type','Papp','category','MDCK','source',source);

                                    %Blood:Plasma Ratio
                                    details.species = 'human';
                                    details.strain  = '';
                                    details.sex     = '';
                                    obj(idx(j)) = obj(idx(j)).addPropertyValue('bp','cmpd',g.CmpdNumber,'value',p.SRBP,'species','human','details',details,'source',sprintf('%s Human',source));

                                    details.species = 'rat';
                                    details.strain = '';
                                    details.sex = '';
                                    obj(idx(j)) = obj(idx(j)).addPropertyValue('bp','cmpd',g.CmpdNumber,'value',p.SRBP_rat,'species','rat','details',details,'source',sprintf('%s Rat',source));

                                    %fup
                                    human = p.SPrUnbnd/100;
                                    obj(idx(j)) = obj(idx(j)).addPropertyValue('fup','cmpd',g.CmpdNumber,'value',human,'species','human','source',sprintf('%s Human',source));

                                    rat = p.SPrUnbnd_rat/100;
                                    obj(idx(j)) = obj(idx(j)).addPropertyValue('fup','cmpd',g.CmpdNumber,'value',rat,'species','rat','source',sprintf('%s Rat',source));
                                    
                                    % fu microsomes
                                    obj(idx(j)) = obj(idx(j)).addPropertyValue('bind','cmpd',g.CmpdNumber,'matrix','microsomes','fu',p.SFumic,'species','human','source',sprintf('%s Human',source));
                                    
                                    %Add CYP/UGT information
                                    obj(idx(j)) = obj(idx(j)).addPropertyValue('enz','Enzyme','CYP 1A2', 'Substrate',m.CYP_1A2_Substr, 'Inhibitor',m.MET_1A2_Inh, 'source',source);
                                    obj(idx(j)) = obj(idx(j)).addPropertyValue('enz','Enzyme','CYP 2A6', 'Substrate',m.CYP_2A6_Substr, 'Inhibitor','N/A',         'source',source);
                                    obj(idx(j)) = obj(idx(j)).addPropertyValue('enz','Enzyme','CYP 2B6', 'Substrate',m.CYP_2B6_Substr, 'Inhibitor','N/A',         'source',source);
                                    obj(idx(j)) = obj(idx(j)).addPropertyValue('enz','Enzyme','CYP 2C8', 'Substrate',m.CYP_2C8_Substr, 'Inhibitor','N/A',         'source',source);
                                    obj(idx(j)) = obj(idx(j)).addPropertyValue('enz','Enzyme','CYP 2C9', 'Substrate',m.CYP_2C9_Substr, 'Inhibitor',m.MET_2C9_Inh, 'source',source);
                                    obj(idx(j)) = obj(idx(j)).addPropertyValue('enz','Enzyme','CYP 2C19','Substrate',m.CYP_2C19_Substr,'Inhibitor',m.MET_2C19_Inh,'source',source);
                                    obj(idx(j)) = obj(idx(j)).addPropertyValue('enz','Enzyme','CYP 2D6', 'Substrate',m.CYP_2D6_Substr, 'Inhibitor',m.MET_2D6_Inh, 'source',source);
                                    obj(idx(j)) = obj(idx(j)).addPropertyValue('enz','Enzyme','CYP 2E1', 'Substrate',m.CYP_2E1_Substr, 'Inhibitor','N/A',         'source',source);
                                    obj(idx(j)) = obj(idx(j)).addPropertyValue('enz','Enzyme','CYP 3A4', 'Substrate',m.CYP_3A4_Substr, 'Inhibitor',m.MET_3A4_Inh, 'source',source);
                                    obj(idx(j)) = obj(idx(j)).addPropertyValue('enz','Enzyme','UGT 1A1', 'Substrate',m.MET_UGT1A1,     'Inhibitor','N/A',         'source',source);
                                    obj(idx(j)) = obj(idx(j)).addPropertyValue('enz','Enzyme','UGT 1A3', 'Substrate',m.MET_UGT1A3,     'Inhibitor','N/A',         'source',source);
                                    obj(idx(j)) = obj(idx(j)).addPropertyValue('enz','Enzyme','UGT 1A4', 'Substrate',m.MET_UGT1A4,     'Inhibitor','N/A',         'source',source);
                                    obj(idx(j)) = obj(idx(j)).addPropertyValue('enz','Enzyme','UGT 1A6', 'Substrate',m.MET_UGT1A6,     'Inhibitor','N/A',         'source',source);
                                    obj(idx(j)) = obj(idx(j)).addPropertyValue('enz','Enzyme','UGT 1A8', 'Substrate',m.MET_UGT1A8,     'Inhibitor','N/A',         'source',source);
                                    obj(idx(j)) = obj(idx(j)).addPropertyValue('enz','Enzyme','UGT 1A9', 'Substrate',m.MET_UGT1A9,     'Inhibitor','N/A',         'source',source);
                                    obj(idx(j)) = obj(idx(j)).addPropertyValue('enz','Enzyme','UGT 1A10','Substrate',m.MET_UGT1A10,    'Inhibitor','N/A',         'source',source);
                                    obj(idx(j)) = obj(idx(j)).addPropertyValue('enz','Enzyme','UGT 2B7', 'Substrate',m.MET_UGT2B7,     'Inhibitor','N/A',         'source',source);
                                    obj(idx(j)) = obj(idx(j)).addPropertyValue('enz','Enzyme','UGT 2B15','Substrate',m.MET_UGT2B15,    'Inhibitor','N/A',         'source',source);
                                    obj(idx(j)) = obj(idx(j)).addPropertyValue('enz','Enzyme','PgP',     'Substrate',p.SPgp_Substr,    'Inhibitor',p.SPgp_Inh,    'source',source);

                                    % Add predicted CLint (microsome)
                                    obj(idx(j)) = obj(idx(j)).addPropertyValue('met','cmpd',g.CmpdNumber,'CL',str2double(m.CYP_HLM_CLint)*(38/1000),'unit','mL/min/g','system','microsomes','species','human','source',sprintf('%s Human unbound Microsomes',source),'details','unbound');
                                    obj(idx(j)) = obj(idx(j)).addPropertyValue('met','cmpd',g.CmpdNumber,'CL',str2double(m.CYP_RLM_CLint)*(38/1000),'unit','mL/min/g','system','microsomes','species','rat','source',sprintf('%s Rat unbound Microsomes',source),'details','unbound');
                                    obj(idx(j)) = obj(idx(j)).addPropertyValue('met','cmpd',g.CmpdNumber,'CL',str2double(m.CYP_HLM_CLint_bound)*(38/1000),'unit','mL/min/g','system','microsomes','species','human','source',sprintf('%s Human bound Microsomes',source),'details','bound');
                                    obj(idx(j)) = obj(idx(j)).addPropertyValue('met','cmpd',g.CmpdNumber,'CL',str2double(m.CYP_RLM_CLint_bound)*(38/1000),'unit','mL/min/g','system','microsomes','species','rat','source',sprintf('%s Rat bound Microsomes',source),'details','bound');
                                    
                                    % Add predicted human Vss
                                    obj(idx(j)) = obj(idx(j)).addPropertyValue('derived_props','type','human Vdss','value',p.SVd,'units','L/kg','source',sprintf('%s Human',source));
                                    
                                    % Add likelihood of brain penetration
                                    obj(idx(j)) = obj(idx(j)).addPropertyValue('derived_props','type','human BBB Filter','value',p.SBBB_Filter,'units','NA','source',sprintf('%s Human',source));
                                    
                                    % Add log of blood brain Kb
                                    obj(idx(j)) = obj(idx(j)).addPropertyValue('derived_props','type','human logKb brain blood','value',p.SLogBB,'units','dimensionless','source',sprintf('%s Human',source));
                                else
                                    obj(idx(j)) = obj.addPropertyValue('fxn_status','service','BuildMasterADMETCompoundProfile','ReturnedValues',0,'comments',admet(i).Comments);
                                end
                            end
                        end
                    end
                end
            end
        end
        
        function obj = getDPSProperties(obj)
            
            dpsCalc = struct.empty;
            
            mdls = {'ClogP_CXN','ClogP_Day','ChromLogP_v1',...              % logP
                    'pKa',...                                               % pKa
                    'PPBRat','PPBHuman',...                                 % plasma protein binding
                    'MDCK2_Perm_pH65','MDCK2_Perm_pH74','PermPeff_3'...     % permeability
                    'Abraham',...                                           % Abraham alpha (H-bond acidity)
                    'CLint_Human','CLint_Rat','CLiv_Rat',...                % clearance
                    'VDSSRat'};                                             % volume of distribution
                
            timer = tic;
            try
                dpsCalc = smtb.webserv.getQSARpropsDPS(strrep(obj.smiles,'#','%23'),obj.gcn);
            catch ME
                obj = obj.addPropertyValue('fxn_status','service','GetQSARPropertiesByCmpd','time_s',timer,'ReturnedValues',0,'comments',ME.message);
            end
            timer = toc(timer);
            
            if isempty(dpsCalc)
                % now try populating MW and MF
                mw = getSelectDPSmodelOutputs(strrep(obj.smiles,'#','%23'),'pis.mw');
                mf = getSelectDPSmodelOutputs(strrep(obj.smiles,'#','%23'),'pis.mf');
                if ~isempty(mw)
                    obj.mw = str2double(mw(strcmpi({mw.PropertyName},'pis.mw')).PropertyValue);
                else
                    obj.mw = NaN;
                end
                if ~isempty(mf)
                    obj.mf = mf(strcmpi({mf.PropertyName},'pis.mf')).PropertyValue;
                else
                    obj.mf = 'Error retrieving MF';
                end
            else
                obj = obj.addPropertyValue('fxn_status','service','GetQSARPropertiesByCmpd','time_s',timer,'FromCache',dpsCalc.FromCache,'ReturnedValues',1,'comments',dpsCalc.Comments);
                obj.mw = dpsCalc.MolWeight;
                obj.mf = dpsCalc.MolFormula;
                
                src = {'QSAR Chemaxon','QSAR Biobyte/Daylight','QSAR ChromLogP v1',...
                       'QSAR Chemaxon',...
                       'QSAR PPB Rat v 1.3','QSAR PPB Human v 1.3',...
                       'QSAR MDCK v1','QSAR MDCK v2','QSAR Peff v3',...
                       'Abraham alpha from SMILES',...
                       'QSAR Met_Intrinsic_Clearance_Human_Microsome_v4','QSAR Rat Intrinsic Clearance Class Microsomes','QSAR Met__in-vivo_clearance_v1.0',...
                       'QSAR Volume_of_Distribution_Rat_v1.2'};
            end
            
            if ~isempty(dpsCalc) && ~contains(dpsCalc.Comments,'No QSAR property')
                % add QSAR ChromLogD model version 9 here. Will update
                % after DPS ADME models are finished
%                 if isnan(str2double(dpsCalc.ChromLogD_v9_cChromlogD_Value))
%                     try
%                         out = getSelectDPSmodelOutputs(strrep(obj.smiles,'#','%23'),'ChromlogD');
%                         obj = obj.addPropertyValue('logP','cmpd',obj.gcn,'value',str2double(out(strcmpi({out.PropertyName},'ChromlogD.value')).PropertyValue),'source','QSAR Chrom LogD pH 7.4 - version 9');
%                     catch
%                     end
%                 else
                if ~isnan(str2double(dpsCalc.ChromLogD_v9_cChromlogD_Value))
                    obj = obj.addPropertyValue('logP','cmpd',obj.gcn,'value',str2double(dpsCalc.ChromLogD_v9_cChromlogD_Value),'source','QSAR Chrom LogD pH 7.4 - version 9');
                end
                if ~isempty(dpsCalc.ChromLogD_v8_dpsErrorString)
                    obj = obj.addPropertyValue('logP','cmpd',obj.gcn,'value',dpsCalc.ChromLogD_v8_CHROM_LOGD_PH74_MEAN,'source','QSAR Chrom LogD pH 7.4 - version 8');
                else
                    obj = obj.addPropertyValue('logP','cmpd',obj.gcn,'value',dpsCalc.ChromlogD_Value,'source','QSAR Chrom LogD pH 7.4 - version 8');
                end
                for i = 1:length(mdls)
                    switch mdls{i}
                        case {'ClogP_CXN','ClogP_Day','ChromLogP_v1','ChromlogD'}
                            if ~isempty(dpsCalc.(sprintf('%s_Value',mdls{i}))) || ~isnan(dpsCalc.(sprintf('%s_Value',mdls{i})))
                                obj = obj.addPropertyValue('logP','cmpd',obj.gcn,'value',dpsCalc.(sprintf('%s_Value',mdls{i})),'source',src{strcmpi(mdls,mdls{i})});
                            else
                                warning('%s DPS model error',src{strcmpi(mdls,mdls{i})})
                            end
                        case 'pKa'
                            if isempty(dpsCalc.(sprintf('%s_ErrorMsg',mdls{i}))) || strcmpi(dpsCalc.(sprintf('%s_ErrorMsg',mdls{i})),'ok')
                                pka_a = [dpsCalc.AcidpKa1 dpsCalc.AcidpKa2 dpsCalc.AcidpKa3 dpsCalc.AcidpKa4];  acidic = pka_a(~(pka_a==14));
                                pka_b = [dpsCalc.BasepKa1 dpsCalc.BasepKa2 dpsCalc.BasepKa3 dpsCalc.BasepKa4];  basic = pka_b(~(pka_b==0));
                                obj = obj.addPropertyValue('pKa','acidic',acidic,'basic',basic,'source',src{strcmpi(mdls,mdls{i})});
                            else
                                warning('%s DPS model error',src{strcmpi(mdls,mdls{i})})
                            end
                        case {'PPBRat','PPBHuman'}
                            if isempty(dpsCalc.(sprintf('%s_ErrorMsg',mdls{i}))) || strcmpi(dpsCalc.(sprintf('%s_ErrorMsg',mdls{i})),'ok')
                                ppb = (100 - dpsCalc.(sprintf('%s_PctB_Value',mdls{i})))/100;
                                obj   = obj.addPropertyValue('fup','cmpd',obj.gcn,'value',ppb,'species',lower(strrep(mdls{i},'PPB','')),'source',src{strcmpi(mdls,mdls{i})});
                            else
                                warning('%s DPS model error',src{strcmpi(mdls,mdls{i})})
                            end
                        case {'MDCK2_Perm_pH65','MDCK2_Perm_pH74'}
                            if isempty(dpsCalc.MDCK2_Perm_ErrorMsg) || strcmpi(dpsCalc.MDCK2_Perm_ErrorMsg,'ok')
                                if contains(mdls{i},'65')
                                    obj = obj.addPropertyValue('perm','cmpd',obj.gcn,'value',dpsCalc.(sprintf('%s_Value',mdls{i})),'unit','nm/s','type','Papp','category','MDCK pH 6.5','source',src{strcmpi(mdls,mdls{i})});
                                else
                                    obj = obj.addPropertyValue('perm','cmpd',obj.gcn,'value',dpsCalc.(sprintf('%s',mdls{i})),'unit','nm/s','type','Papp','category','MDCK pH 7.4','source',src{strcmpi(mdls,mdls{i})});
                                end
                            else
                                warning('%s DPS model error',src{strcmpi(mdls,mdls{i})})
                            end
                        case 'Abraham'
                            if isempty(dpsCalc.(sprintf('%s_ErrorMsg',mdls{i}))) || strcmpi(dpsCalc.(sprintf('%s_ErrorMsg',mdls{i})),'ok')
                                obj   = obj.addPropertyValue('descriptors','cmpd',obj.gcn,'value',str2double(dpsCalc.(sprintf('%s_Alpha',mdls{i}))),'type','HbondA','source','Abraham alpha from SMILES');
                            else
                                warning('%s DPS model error',src{strcmpi(mdls,mdls{i})})
                            end
                        case {'CLint_Human','CLint_Rat'}
                            if isempty(dpsCalc.(sprintf('%s_ErrorMsg',mdls{i}))) || strcmpi(dpsCalc.(sprintf('%s_ErrorMsg',mdls{i})),'ok')
                                sp = lower(strsplit(mdls{i},'_')); sp = sp{2};
                                if strcmpi(sp,'human')
                                    obj   = obj.addPropertyValue('met','cmpd',obj.gcn,'CL',dpsCalc.(sprintf('%s_Value',mdls{i})),'unit','mL/min/g','system','microsomes','species',sp,'source',src{strcmpi(mdls,mdls{i})},'details','bound');
                                    obj   = obj.addPropertyValue('derived_props','type','HLM clearance class','value',dpsCalc.(sprintf('%s_Class',mdls{i})),'units','mL/min/g','source',src{strcmpi(mdls,mdls{i})});
                                else
                                    obj   = obj.addPropertyValue('derived_props','type','rat intrinsic clearance class','value',dpsCalc.(sprintf('%s_Class',mdls{i})),'units','mL/min/g','source',src{strcmpi(mdls,mdls{i})});
                                    obj   = obj.addPropertyValue('met','cmpd',obj.gcn,'CL',dpsCalc.(sprintf('%s_Value',mdls{i})),'unit','mL/min/g','system','microsomes','species','rat','source',src{strcmpi(mdls,mdls{i})},'details','bound');
                                end 
                            else
                                warning('%s DPS model error',src{strcmpi(mdls,mdls{i})})
                            end
                        case 'CLiv_Rat'
                            if isempty(dpsCalc.(sprintf('%s_ErrorMsg',mdls{i}))) || strcmpi(dpsCalc.(sprintf('%s_ErrorMsg',mdls{i})),'ok')
                                obj   = obj.addPropertyValue('derived_props','type','rat in vivo clearance class','value',dpsCalc.(sprintf('%s_Class',mdls{i})),'units','mL/min/kg','source',src{strcmpi(mdls,mdls{i})});
                                obj   = obj.addPropertyValue('derived_props','type','rat in vivo clearance','value',10^dpsCalc.(sprintf('%s_LogCL_Score',mdls{i})),'units','mL/min/kg','source',src{strcmpi(mdls,mdls{i})});
                            else
                                warning('%s DPS model error',src{strcmpi(mdls,mdls{i})})
                            end
                        case 'VDSSRat'
                            if isempty(dpsCalc.(sprintf('%s_ErrorMsg',mdls{i}))) || strcmpi(dpsCalc.(sprintf('%s_ErrorMsg',mdls{i})),'ok')
                                obj   = obj.addPropertyValue('derived_props','type','rat Vdss','value',dpsCalc.(sprintf('%s_Value',mdls{i})),'units','L/kg','source',src{strcmpi(mdls,mdls{i})});
                            else
                                warning('%s DPS model error',src{strcmpi(mdls,mdls{i})})
                            end
                    end
                end
            else
                warning('Error retrieving DPS QSAR properties.')
            end
        end
        
        function obj = getAdamantisProperties(obj,mtbc)
            
            prms = ['ChromLogP_v1.value;clogp_cxn.value;clogp_day.value;ChromlogD.value;'... %logP
                    'pKacxn.a1;pKacxn.a2;pKacxn.a3;pKacxn.a4;'...            %pKa - acidic
                    'pKacxn.b1;pKacxn.b2;pKacxn.b3;pKacxn.b4;'...            %pKa - basic
                    'PPBrat_v1_3.PctB;PPBhuman_v1_3.PctB;'...                %Plasma Protein Binding
                    'MDCK2_pH74_v2.Perm_nm_s;MDCK2.Perm_pH65_nm_sec;permeability_v3.Peff;'... %Permeability
                    'abraham.alpha;'...                                      %Abrahams alpha (H-bond acidity)
                    'clint_v4.Value_ml_min_g;clint_v4.Class;Clr_Rat.Clr_Rat_Class;Clr_Rat.Clr_Rat_value;CLiv_rat_v1.CL_class;CLiv_rat_v1.logCL_score;'...       %clearance
                    'VDSSrat_v1_2.VDSS;'...                                  %Vss
                    'ChromLogP_v1.dpsErrorString;CXN.dpsErrorString;clogp_day.dpsErrorString;ChromlogD.dpsErrorString;CXN_pKa.dpsErrorString;PPBrat_v1_3.dpsErrorString;PPBhuman_v1_3.dpsErrorString;MDCK2.dpsErrorString;permeability_v3.dpsErrorString;abraham.alpha.dpsErrorString;clint_v4.dpsErrorString;Clr_Rat.dpsErrorString;CLiv_rat_v1.dpsErrorString;VDSSrat_v1_2.dpsErrorString;'... % check for model errors
                    'ChromlogD.version']; % add model versions where appropriate
            dpsCalc = GetSelectCalculatedDerivedPropertiesDPS(mtbc, prms, obj.smiles);
            props = dpsCalc.CmpdProperties;
            
            warning OFF BACKTRACE
            %logP
            if isempty(props(27).PropertyValue) || strcmpi(props(27).PropertyValue,'ok')
                obj = obj.addPropertyValue('logP','cmpd',obj.gcn,'value',str2double(props(2).PropertyValue),'source','QSAR ChromLogP v1');
            else
                warning('ChromLogP model error : "%s"',props(27).PropertyValue)
            end
            
            if isempty(props(28).PropertyValue) || strcmpi(props(28).PropertyValue,'ok')
                obj = obj.addPropertyValue('logP','cmpd',obj.gcn,'value',str2double(props(3).PropertyValue),'source','QSAR Chemaxon');
            else
                warning('ChemAxon logP model error : "%s"',props(28).PropertyValue)
            end
            
            if isempty(props(29).PropertyValue) || strcmpi(props(29).PropertyValue,'ok')
                obj = obj.addPropertyValue('logP','cmpd',obj.gcn,'value',str2double(props(4).PropertyValue),'source','QSAR Biobyte/Daylight');
            else
                warning('Bioyte/Daylight logP model error : "%s"',props(29).PropertyValue)
            end
            
            if isempty(props(30).PropertyValue) || strcmpi(props(30).PropertyValue,'ok')
                obj = obj.addPropertyValue('logP','cmpd',obj.gcn,'value',str2double(props(5).PropertyValue),'source',sprintf('QSAR Chrom LogD pH 7.4 - version %s',props(41).PropertyValue));
            else
                warning('Chrom logD model error : "%s"',props(30).PropertyValue)
            end
            
            %pKa
            if isempty(props(31).PropertyValue) || strcmpi(props(31).PropertyValue,'ok')
                pka_a = str2double({props(6:9).PropertyValue});  acidic = pka_a(~(pka_a==14));
                pka_b = str2double({props(10:13).PropertyValue}); basic = pka_b(~(pka_b==0));
                obj = obj.addPropertyValue('pKa','acidic',acidic,'basic',basic,'source','QSAR Chemaxon');
            else
                warning('ChemAxon pKa model error : "%s"',props(31).PropertyValue)
            end
            
            %fup
            if isempty(props(32).PropertyValue) || strcmpi(props(32).PropertyValue,'ok')
                rat = (100-str2double(props(14).PropertyValue))/100;
                obj = obj.addPropertyValue('fup','cmpd',obj.gcn,'value',rat,'species','rat','source','QSAR PPB Rat v 1.3');
            else
                warning('PPB Rat v1.3 model error : "%s"',props(32).PropertyValue)
            end
            
            if isempty(props(33).PropertyValue) || strcmpi(props(33).PropertyValue,'ok')
                human = (100-str2double(props(15).PropertyValue))/100;
                obj   = obj.addPropertyValue('fup','cmpd',obj.gcn,'value',human,'species','human','source','QSAR PPB Human v 1.3');
            else
                warning('PPB Human v1.3 model error : "%s"',props(33).PropertyValue)
            end 
            
            %Permeability
            if isempty(props(34).PropertyValue) || strcmpi(props(34).PropertyValue,'ok')
                obj = obj.addPropertyValue('perm','cmpd',obj.gcn,'value',str2double(props(16).PropertyValue),'unit','nm/s','type','MDCK pH 7.4','source','QSAR MDCK v2');
                obj = obj.addPropertyValue('perm','cmpd',obj.gcn,'value',str2double(props(17).PropertyValue),'unit','nm/s','type','MDCK pH 6.5','source','QSAR MDCK v1');
            else
                warning('MDCK2 model error: "%s"',props(34).PropertyValue)
            end
            
            if isempty(props(35).PropertyValue) || strcmpi(props(35).PropertyValue,'ok')
                obj = obj.addPropertyValue('perm','cmpd',obj.gcn,'value',str2double(props(18).PropertyValue),'unit','nm/s','type','Peff','source','QSAR Peff v3');
            else
                warning('Peff v3 model error: "%s"',props(35).PropertyValue)
            end
            
            %Abrahams alpha
            if isempty(props(36).PropertyValue) || strcmpi(props(36).PropertyValue,'ok')
                obj   = obj.addPropertyValue('descriptors','cmpd',obj.gcn,'value',str2double(props(19).PropertyValue),'type','HbondA','source','Abraham alpha from SMILES');
            else
                warning('Abraham alpha error : "%s"',props(36).PropertyValue)
            end 
            
            %clearance
            if isempty(props(37).PropertyValue) || strcmpi(props(37).PropertyValue,'ok')
                obj   = obj.addPropertyValue('met','cmpd',obj.gcn,'CL',str2double(props(20).PropertyValue),'unit','mL/min/g','system','microsomes','species','human','source','QSAR Met_Intrinsic_Clearance_Human_Microsome_v4','details','bound');
                obj   = obj.addPropertyValue('derived_props','type','HLM clearance class','value',props(21).PropertyValue,'units','mL/min/g','source','QSAR Met_Intrinsic_Clearance_Human_Microsome_v4');
            else
                warning('clint_v4 model error : "%s"',props(37).PropertyValue)
            end

            if isempty(props(38).PropertyValue) || strcmpi(props(38).PropertyValue,'ok')
                obj   = obj.addPropertyValue('derived_props','type','rat intrinsic clearance class','value',props(22).PropertyValue,'units','mL/min/g','source','QSAR Rat Intrinsic Clearance Class');
                obj   = obj.addPropertyValue('met','cmpd',obj.gcn,'CL',str2double(props(23).PropertyValue),'unit','mL/min/g','system','microsomes','species','rat','source','QSAR Rat Intrinsic Clearance Class Microsomes','details','bound');
            else
                warning('Clr_Rat model error : "%s"',props(38).PropertyValue)
            end
            
            if isempty(props(39).PropertyValue) || strcmpi(props(39).PropertyValue,'ok')
                obj   = obj.addPropertyValue('derived_props','type','rat in vivo clearance class','value',props(24).PropertyValue,'units','mL/min/kg','source','QSAR Met__in-vivo_clearance_v1.0');
                obj   = obj.addPropertyValue('derived_props','type','rat in vivo clearance','value',10^str2double(props(25).PropertyValue),'units','mL/min/kg','source','QSAR Met__in-vivo_clearance_v1.0');
            else
                warning('CLiv_rat_v1 model error : "%s"',props(39).PropertyValue)
            end
            
            % Vss
            if isempty(props(40).PropertyValue) || strcmpi(props(40).PropertyValue,'ok')
                obj   = obj.addPropertyValue('derived_props','type','rat Vdss','value',str2double(props(26).PropertyValue),'units','L/kg','source','QSAR Volume_of_Distribution_Rat_v1.2');
            else
                warning('VDSSrat_v1_2 : "%s"',props(40).PropertyValue)
            end
            
            
            warning ON BACKTRACE
            
        end
        
        function obj = getKATEProperties(obj)
            timer = tic;
            p = getKATEDataJSON(obj.name);
            timer = toc(timer);
            obj = obj.addPropertyValue('fxn_status','service','getKATEDataJSON','time_s',timer);
            if isempty(p.PhysChemSumLST)  % if a generic name is input, p entries will be empty. this runs it with gcn if that is the case.
                p = getKATEDataJSON(obj.gcn);
            end
            if ~isempty(p.PhysChemSumLST)
                [~,idx] = sort({p.PhysChemSumLST.CmpdNumber});
                p.PhysChemSumLST = p.PhysChemSumLST(idx);
                for i = 1:length(p.PhysChemSumLST)
                    pchem = p.PhysChemSumLST(i);
                    %LogP
                    if ~isempty(pchem.CHROM_LOGP),           obj = obj.addPropertyValue('logP','value',str2double(pchem.CHROM_LOGP),'source','KATE Chrom logP','cmpd',pchem.CmpdNumber); end
                    if ~isempty(pchem.CHI_LOGP_MEAN),        obj = obj.addPropertyValue('logP','value',str2double(pchem.CHI_LOGP_MEAN),'source','KATE CHI logP','cmpd',pchem.CmpdNumber); end
                    if ~isempty(pchem.LOGK_IAM),             obj = obj.addPropertyValue('logP','value',0.29.*exp(str2double(pchem.LOGK_IAM))+0.70,  'source','KATE LogK IAM, adjusted to logP','cmpd',pchem.CmpdNumber); end
                    if ~isempty(pchem.CHROM_LOGD_PH74_MEAN), obj = obj.addPropertyValue('logP','value',str2double(pchem.CHROM_LOGD_PH74_MEAN),'source','KATE Chrom LogD pH 7.4','cmpd',pchem.CmpdNumber); end
                    if ~isempty(pchem.CHI_LOGD_PH74_MEAN),   obj = obj.addPropertyValue('logP','value',str2double(pchem.CHI_LOGD_PH74_MEAN),'source','KATE CHI LogD pH 7.4','cmpd',pchem.CmpdNumber); end
                    if ~isempty(pchem.CHI_IAM_MEAN),         obj = obj.addPropertyValue('logP','value',str2double(pchem.CHI_IAM_MEAN),'source','KATE CHI IAM','cmpd',pchem.CmpdNumber); end
                    %Permeability
                    if ~isempty(pchem.PERM_NUM_PH705_MEAN), obj = obj.addPropertyValue('perm','Modifier',pchem.PERM_NUM_MOD_PH705,'value',str2double(pchem.PERM_NUM_PH705_MEAN),'unit','nm/s','type','Ppassive','category','AMP - pH 7.05','source','KATE','cmpd',pchem.CmpdNumber); end
                    if ~isempty(pchem.PERM_NUM_PH7_4_MEAN), obj = obj.addPropertyValue('perm','Modifier',pchem.PERM_NUM_MOD_PH7_4,'value',str2double(pchem.PERM_NUM_PH7_4_MEAN),'unit','nm/s','type','Ppassive','category','AMP - pH 7.4','source','KATE','cmpd',pchem.CmpdNumber); end
                    [obj.perm(strcmp({obj.perm.Modifier},'=')).Modifier] = deal([]);
                end
            elseif ~isempty(p.CHI_SLST) %This is essentially a catch in case there is an issue with the PHYSCHEM_SUMMARY table.  Issue has been resolved with few compounds that have been tested but will keep code.
                [~,idx] = sort({p.CHI_SLST.CmpdNumber});
                p.CHI_SLST = p.CHI_SLST(idx);
                for i = 1:length(p.CHI_SLST)
                    pchem = p.CHI_SLST(i);
                    if ~isempty(pchem.AVG_CHROM_LOGP),          obj = obj.addPropertyValue('logP','value',str2double(pchem.AVG_CHROM_LOGP),'source','KATE Chrom logP','cmpd',pchem.CmpdNumber); end
                    if ~isempty(pchem.AVG_CHI_LOGP),            obj = obj.addPropertyValue('logP','value',str2double(pchem.AVG_CHI_LOGP),'source','KATE CHI logP','cmpd',pchem.CmpdNumber); end
                    if ~isempty(pchem.AVG_CHROM_LOGD_PH_7_4),   obj = obj.addPropertyValue('logP','value',str2double(pchem.AVG_CHROM_LOGD_PH_7_4),'source','KATE Chrom LogD pH 7.4','cmpd',pchem.CmpdNumber); end
                    if ~isempty(pchem.AVG_CHI_LOGD_PH_7_4),     obj = obj.addPropertyValue('logP','value',str2double(pchem.AVG_CHI_LOGD_PH_7_4),'source','KATE CHI LogD pH 7.4','cmpd',pchem.CmpdNumber); end
                    if ~isempty(pchem.AVG_CHI_IAM),             obj = obj.addPropertyValue('logP','value',str2double(pchem.AVG_CHI_IAM),'source','KATE CHI IAM','cmpd',pchem.CmpdNumber); end
                end
            end
            
            %CHIN LogP
            if ~isempty(p.CHI_SLST)
                [~,idx] = sort({p.CHI_SLST.CmpdNumber});
                p.CHI_SLST = p.CHI_SLST(idx);
                for i = 1:length(p.CHI_SLST)
                   chi_s = p.CHI_SLST(i);
                   if ~isempty([obj.descriptors(strcmpi({obj.descriptors.type},'HbondA')).value])
                       CHIN = max([str2double(chi_s.AVG_CHI_PH_10_5) str2double(chi_s.AVG_CHI_PH_2_0) str2double(chi_s.AVG_CHI_PH_7_4)]);
                       A = obj.descriptors(strcmpi({obj.descriptors.type},'HbondA')).value;
                       CHINlogP = 0.054*CHIN + 1.319*A - 1.877;
                       obj = obj.addPropertyValue('logP','value',CHINlogP,'source','KATE CHIN logP','cmpd',chi_s.CmpdNumber);
                   end
                end
            end
            
            %Permeability
            if ~isempty(p.InVitroPermLST)
                [~,idx] = sort({p.InVitroPermLST.CmpdNumber});
                p.InVitroPermLST = p.InVitroPermLST(idx);
                for i = 1:length(p.InVitroPermLST)
                    perm = p.InVitroPermLST(i);
                    if ~isempty(perm.CMPD_PAPP_NM_SEC),         obj = obj.addPropertyValue('perm','Modifier',perm.CMPD_PAPP_MOD,'value',str2double(perm.CMPD_PAPP_NM_SEC),'efflux_ratio',str2double(perm.EFFLUX_RATIO),'pct_recovery',str2double(perm.PCT_RECOVERY),'unit','nm/s','type',sprintf('%s %s pH %s/%s',perm.CellLine,perm.MEAS_DIRECTION,perm.DONOR_PH,perm.RECEIVER_PH),'source','KATE','cmpd',perm.CmpdNumber,'comments',perm.Comments,'study',perm.StudyNumber); end
                    if ~isempty(perm.CMPD_PAPP_INH_NM_SEC),     obj = obj.addPropertyValue('perm','Modifier',perm.CMPD_PAPP_INH_MOD,'value',str2double(perm.CMPD_PAPP_INH_NM_SEC),'efflux_ratio',str2double(perm.EFFLUX_RATIO_INH),'pct_recovery',str2double(perm.PCT_RECOVERY_INH),'unit','nm/s','type',sprintf('%s %s pH %s/%s Inh %s',perm.CellLine,perm.MEAS_DIRECTION,perm.DONOR_PH,perm.RECEIVER_PH,perm.INHIBITOR),'source','KATE','cmpd',perm.CmpdNumber,'comments',perm.Comments,'study',perm.StudyNumber); end
                end
            end
            
            %LogD
            
            %pKa
            
            %Protein Binding
            if ~isempty(p.ProteinBindingLST)
                bnd = struct('cmpd',    {p.ProteinBindingLST.CmpdNumber},...
                             'species', [],...
                             'strain',  [],...
                             'sex',     [],...
                             'matrix',  {p.ProteinBindingLST.MATRIX},...
                             'method',  {p.ProteinBindingLST.METHOD},...
                             'fu_mod',  {p.ProteinBindingLST.CMPD_PCT_BINDING_MOD},...
                             'fu_pct',  num2cell(str2double({p.ProteinBindingLST.FU_PERCENT})./100),...
                             'conc',    num2cell(str2double({p.ProteinBindingLST.FU_PERCENT_MOD})),...
                             'conc_units',{p.ProteinBindingLST.CMPD_CONC_UNITS},...
                             'comments',{p.ProteinBindingLST.Comments});
                for i = 1:length(p.ProteinBindingLST)
                    bnd(i).species = lower(p.ProteinBindingLST(i).PBanimal.Species);
                    bnd(i).strain  = lower(p.ProteinBindingLST(i).PBanimal.Strain);
                    bnd(i).sex     = lower(p.ProteinBindingLST(i).PBanimal.Gender);
                end
                
                sp = unique(lower({bnd.species}));
                ma = unique(lower({bnd.matrix}));
                cn = unique({bnd.cmpd});
                for i = 1:length(sp)
                    for j = 1:length(ma)
                        for k = 1:length(cn)
                            details = bnd(strcmpi({bnd.species},sp{i}) & strcmpi({bnd.matrix},ma{j}) & strcmpi({bnd.cmpd},cn{k}));
                            if ~isempty(details)
                                for l = 1:length(details)
                                    obj = obj.addPropertyValue('bind','cmpd',cn{k},'matrix',ma{j},'fu',details(l).fu_pct,'species',sp{i},'details',details(l),'source','KATE','method',details(l).method,'comments',details(l).comments);
                                end
                            end
                        end
                    end
                end
            end
            
            %fup
            bnd = obj.bind(strcmpi({obj.bind.source},'kate'));
            sp = unique({bnd.species});
            cn = unique({bnd.cmpd});
            for i = 1:length(sp)
                for j = 1:length(cn)
                    idx = find(strcmp({bnd.species},sp{i}) & strcmp({bnd.cmpd},cn{j}) & strcmpi({bnd.matrix},'plasma'));
                    if ~isempty(idx)
                        for k = 1:length(bnd(idx))
                            obj = obj.addPropertyValue('fup','cmpd',cn{j},'value',bnd(idx(k)).fu,'species',sp{i},'details',bnd(idx(k)).details,'source','KATE','method',bnd(idx(k)).method,'comments',bnd(idx(k)).comments);
                        end
                    end
                end
            end
            
            
            %B:P ratio
            if ~isempty(p.BloodPlasmaPartitionLST)
                b_p = p.BloodPlasmaPartitionLST;
                for i = 1:length(b_p)
                    details = struct('Species',[],'Strain',[],'Sex',[],'Concentration',[],'Units',[],'BP_Ratio',[]);
                    details(1).Concentration = str2double(b_p(i).CONC1_NG_ML);   details(1).BP_Ratio = str2double(b_p(i).B_P_RATIO1);
                    details(2).Concentration = str2double(b_p(i).CONC2_NG_ML);   details(2).BP_Ratio = str2double(b_p(i).B_P_RATIO2);
                    details(3).Concentration = str2double(b_p(i).CONC3_NG_ML);   details(3).BP_Ratio = str2double(b_p(i).B_P_RATIO3);
                    [details.species] = deal(lower([b_p(i).BPPAnimal.Species]));
                    [details.strain]  = deal(lower([b_p(i).BPPAnimal.Strain]));
                    [details.sex]     = deal(b_p(i).BPPAnimal.Gender);
                    if std([details.BP_Ratio]) > (0.25 * mean([details.BP_Ratio])), flag = 'Standard deviation of B:P Ratio greater than 25% of the mean. This may indicate concentration dependence'; else, flag = ''; end
                    obj = obj.addPropertyValue('bp','cmpd',b_p(i).CmpdNumber,'value',nanmean([details.BP_Ratio]),'species',lower(b_p(i).BPPAnimal.Species),'details',details,'source','KATE','flag',flag,'comments',b_p(i).Comments);
                end
            end
            
            %Intrinsic clearance
            if ~isempty(p.InVitroMetStabLST)
                ms = p.InVitroMetStabLST;
                for i = 1:length(ms)
                    details = struct('species',lower(ms(i).MetStabilityAnimal.Species),'strain',lower(ms(i).MetStabilityAnimal.Strain),'sex',lower(ms(i).MetStabilityAnimal.Gender),'CmpdConc',ms(i).COMD_CONC_UM,'CmpdConcUnit','uM','SystemConc',ms(i).INVITRO_SYSTEM_CONC,'SystemUnit',ms(i).INVITRO_SYSTEM_UNIT);
                    obj = obj.addPropertyValue('met','cmpd',ms(i).CmpdNumber,'Modifier',ms(i).CL_ML_MIN_G_TISSUE_MOD,'CL',str2double(ms(i).CL_ML_MIN_G_TISSUE),'unit','mL/min/g','species',lower(ms(i).MetStabilityAnimal.Species),'system',lower(ms(i).IN_VITRO_SYSTEM),'details',details,'source','KATE');
                end
            end
            
            %PK Studies
            if ~isempty(p.InVivoPharmacokineticsLST)
                pk = p.InVivoPharmacokineticsLST;
                for i = 1:length(pk)
                    obj = obj.addPropertyValue('pk_kate','cmpd',            pk(i).CmpdNumber,...
                                                         'study',           pk(i).StudyNumber,...
                                                         'study_design',    lower(pk(i).STUDY_DESIGN),...
                                                         'species',         lower(pk(i).PKanimal.Species),...
                                                         'strain',          lower(pk(i).PKanimal.Strain),...
                                                         'sex',             lower(pk(i).PKanimal.Gender),...
                                                         'n',               str2double(pk(i).TOT_NBR_ANIMALS),...
                                                         'dose_mgkg',       str2double(pk(i).DOSE_MG_KG),...
                                                         'route',           pk(i).ROUTE_ADMIN,...
                                                         'formulation',     pk(i).FORMULATION,...
                                                         'matrix',          lower(pk(i).SAMPLE_MATRIX),...
                                                         'F',               str2double(pk(i).F_PCT),...
                                                         'F_std',           str2double(pk(i).F_PCT_SD),...
                                                         'Cmax_ngml',       str2double(pk(i).CMAX_NG_ML),...
                                                         'Cmax_sd',         str2double(pk(i).CMAX_SD),...
                                                         'Tmax_hr',         str2double(pk(i).TMAX_HR),...
                                                         'Tmax_sd',         str2double(pk(i).TMAX_SD),...
                                                         'AUC_Inf_nghrmL',  str2double(pk(i).AUC_0INF_NG_HR_ML),...
                                                         'AUC_Inf_std',     str2double(pk(i).AUC_0INF_SD),...
                                                         'TBdy_CL_mlmin_kg',str2double(pk(i).TOT_BDY_CLR_ML_MIN_KG),...
                                                         'TBdy_CL_std',     str2double(pk(i).TOT_BDY_CLR_SD),...
                                                         'Vdss_Lkg',        str2double(pk(i).VDSS_L_KG),...
                                                         'Vdss_std',        str2double(pk(i).VDSS_SD),...
                                                         'T1_2_hr',         str2double(pk(i).T1_2_HR),...
                                                         'T1_2_std',        str2double(pk(i).T1_2_SD),...
                                                         'project',         pk(i).ProjectName,...
                                                         'pk_start_date',   pk(i).ExperimentDate,...
                                                         'pk_date_record',  pk(i).DateRecordCreated,...
                                                         'pk_date_updated', pk(i).DateRecordUpdated,...
                                                         'comments',        pk(i).Comments,...
                                                         'source',          'KATE');
                    flds = fields(obj.pk_kate(end));                
                    for j = 1:length(flds)
                        if prod(isnan(obj.pk_kate(i).(flds{j}))), obj.pk_kate(i).(flds{j}) = []; end
                    end
                end
                [~,idx] = sort([obj.pk_kate.dose_mgkg]);obj.pk_kate = obj.pk_kate(idx);
                [~,idx] = sort({obj.pk_kate.species});  obj.pk_kate = obj.pk_kate(idx);
                [~,idx] = sort({obj.pk_kate.cmpd});     obj.pk_kate = obj.pk_kate(idx);
            end 
        end
        
        function obj = getKATEinVitroProperties(obj)
            url 	= 'https://ivivtapps-prod.gsk.com/KATESearch/KATESearch.svc/REST/';
            format  = 'json';
            opt     = weboptions('Timeout',30,'CertificateFilename','');
            
            timer = tic;
                PhysChemSumLST = webread([url format '/' 'GetPhysChemSummaryResults' '?key=0&cmpd=' obj.name],opt);
                if isempty(PhysChemSumLST), PhysChemSumLST = webread([url format '/' 'GetPhysChemSummaryResults' '?key=0&cmpd=' obj.gcn],opt); end

                CHI_SLST = webread([url format '/' 'GetCHI_S' '?key=0&cmpd=' obj.name],opt);
                if isempty(CHI_SLST), CHI_SLST = webread([url format '/' 'GetCHI_S' '?key=0&cmpd=' obj.gcn],opt); end

                InVitroPermLST = webread([url format '/' 'GetInVitroPermResults' '?key=0&cmpd=' obj.name],opt);
                if isempty(InVitroPermLST), InVitroPermLST = webread([url format '/' 'GetInVitroPermResults' '?key=0&cmpd=' obj.gcn],opt); end
                
                ProteinBindingLST = webread([url format '/' 'GetProteinBindingResults' '?key=0&cmpd=' obj.name],opt);
                if isempty(ProteinBindingLST), ProteinBindingLST = webread([url format '/' 'GetProteinBindingResults' '?key=0&cmpd=' obj.gcn],opt); end
                
                BloodPlasmaPartitionLST = webread([url format '/' 'GetBloodPlasmaPartitioningResults' '?key=0&cmpd=' obj.name],opt);
                if isempty(BloodPlasmaPartitionLST), BloodPlasmaPartitionLST = webread([url format '/' 'GetBloodPlasmaPartitioningResults' '?key=0&cmpd=' obj.gcn],opt); end
                
                InVitroMetStabLST = webread([url format '/' 'GetInVitroMetStabResults' '?key=0&cmpd=' obj.name],opt);
                if isempty(InVitroMetStabLST), InVitroMetStabLST = webread([url format '/' 'GetInVitroMetStabResults' '?key=0&cmpd=' obj.gcn],opt); end
                
                pka = smtb.webserv.getKATEpKa(obj.name);
                if isempty(pka), pka = smtb.webserv.getKATEpKa(obj.gcn); end
                
                sol = smtb.webserv.getKATEsolubility(obj.name);
                if isempty(sol), sol = smtb.webserv.getKATEsolubility(obj.gcn); end
                
            timer = toc(timer);
            if isempty(PhysChemSumLST) && isempty(CHI_SLST) && isempty(InVitroPermLST) && isempty(ProteinBindingLST) && isempty(BloodPlasmaPartitionLST) && isempty(InVitroMetStabLST) && isempty(pka) && isempty(sol)
                obj = obj.addPropertyValue('fxn_status','service','getKATEinVitroProperties','time_s',timer,'ReturnedValues',0,'comments',sprintf('No in vitro data found for %s',obj.name));
            else
                obj = obj.addPropertyValue('fxn_status','service','getKATEinVitroProperties','time_s',timer,'ReturnedValues',1);
            end
            if ~isempty(PhysChemSumLST)
                [~,idx] = sort({PhysChemSumLST.CmpdNumber});
                PhysChemSumLST = PhysChemSumLST(idx);
                for i = 1:length(PhysChemSumLST)
                    pchem = PhysChemSumLST(i);
                    %LogP
                    if ~isempty(pchem.CHROM_LOGP),           obj = obj.addPropertyValue('logP','value',str2double(pchem.CHROM_LOGP),'source','KATE Chrom logP','cmpd',pchem.CmpdNumber,'lnb_ref',pchem.ExperimentLNBRef,'prep_lnb_ref',pchem.PrepLNBRef,'date_experiment',pchem.ExperimentDate,'date_record',pchem.ExperimentDate,'date_updated',pchem.DateRecordUpdated); end
                    if ~isempty(pchem.CHI_LOGP_MEAN),        obj = obj.addPropertyValue('logP','value',str2double(pchem.CHI_LOGP_MEAN),'source','KATE CHI logP','cmpd',pchem.CmpdNumber,'lnb_ref',pchem.ExperimentLNBRef,'prep_lnb_ref',pchem.PrepLNBRef,'date_experiment',pchem.ExperimentDate,'date_record',pchem.ExperimentDate,'date_updated',pchem.DateRecordUpdated); end
                    if ~isempty(pchem.LOGK_IAM),             obj = obj.addPropertyValue('logP','value',0.29.*exp(str2double(pchem.LOGK_IAM))+0.70,  'source','KATE LogK IAM, adjusted to logP','cmpd',pchem.CmpdNumber,'lnb_ref',pchem.ExperimentLNBRef,'prep_lnb_ref',pchem.PrepLNBRef,'date_experiment',pchem.ExperimentDate,'date_record',pchem.ExperimentDate,'date_updated',pchem.DateRecordUpdated); end
                    if ~isempty(pchem.CHROM_LOGD_PH74_MEAN), obj = obj.addPropertyValue('logP','value',str2double(pchem.CHROM_LOGD_PH74_MEAN),'source','KATE Chrom LogD pH 7.4','cmpd',pchem.CmpdNumber,'lnb_ref',pchem.ExperimentLNBRef,'prep_lnb_ref',pchem.PrepLNBRef,'date_experiment',pchem.ExperimentDate,'date_record',pchem.ExperimentDate,'date_updated',pchem.DateRecordUpdated); end
                    if ~isempty(pchem.CHI_LOGD_PH74_MEAN),   obj = obj.addPropertyValue('logP','value',str2double(pchem.CHI_LOGD_PH74_MEAN),'source','KATE CHI LogD pH 7.4','cmpd',pchem.CmpdNumber,'lnb_ref',pchem.ExperimentLNBRef,'prep_lnb_ref',pchem.PrepLNBRef,'date_experiment',pchem.ExperimentDate,'date_record',pchem.ExperimentDate,'date_updated',pchem.DateRecordUpdated); end
                    if ~isempty(pchem.CHI_IAM_MEAN),         obj = obj.addPropertyValue('logP','value',str2double(pchem.CHI_IAM_MEAN),'source','KATE CHI IAM','cmpd',pchem.CmpdNumber,'lnb_ref',pchem.ExperimentLNBRef,'prep_lnb_ref',pchem.PrepLNBRef,'date_experiment',pchem.ExperimentDate,'date_record',pchem.ExperimentDate,'date_updated',pchem.DateRecordUpdated); end
                    %Permeability
                    if ~isempty(pchem.PERM_NUM_PH705_MEAN), obj = obj.addPropertyValue('perm','Modifier',pchem.PERM_NUM_MOD_PH705,'value',str2double(pchem.PERM_NUM_PH705_MEAN),'unit','nm/s','type','Ppassive','category','AMP - pH 7.05','source','KATE','cmpd',pchem.CmpdNumber,'lnb_ref',pchem.ExperimentLNBRef,'prep_lnb_ref',pchem.PrepLNBRef,'date_experiment',pchem.ExperimentDate,'date_record',pchem.ExperimentDate,'date_updated',pchem.DateRecordUpdated); end
                    if ~isempty(pchem.PERM_NUM_PH7_4_MEAN), obj = obj.addPropertyValue('perm','Modifier',pchem.PERM_NUM_MOD_PH7_4,'value',str2double(pchem.PERM_NUM_PH7_4_MEAN),'unit','nm/s','type','Ppassive','category','AMP - pH 7.4','source','KATE','cmpd',pchem.CmpdNumber,'lnb_ref',pchem.ExperimentLNBRef,'prep_lnb_ref',pchem.PrepLNBRef,'date_experiment',pchem.ExperimentDate,'date_record',pchem.ExperimentDate,'date_updated',pchem.DateRecordUpdated); end
                    [obj.perm(strcmp({obj.perm.Modifier},'=')).Modifier] = deal([]);
                end
            elseif ~isempty(CHI_SLST) %This is essentially a catch in case there is an issue with the PHYSCHEM_SUMMARY table.  Issue has been resolved with few compounds that have been tested but will keep code.
                [~,idx] = sort({CHI_SLST.CmpdNumber});
                CHI_SLST = CHI_SLST(idx);
                for i = 1:length(CHI_SLST)
                    pchem = CHI_SLST(i);
                    if ~isempty(pchem.AVG_CHROM_LOGP),          obj = obj.addPropertyValue('logP','value',str2double(pchem.AVG_CHROM_LOGP),'source','KATE Chrom logP','cmpd',pchem.CmpdNumber,'lnb_ref',pchem.ExperimentLNBRef,'prep_lnb_ref',pchem.PrepLNBRef,'date_experiment',pchem.ExperimentDate,'date_record',pchem.ExperimentDate,'date_updated',pchem.DateRecordUpdated); end
                    if ~isempty(pchem.AVG_CHI_LOGP),            obj = obj.addPropertyValue('logP','value',str2double(pchem.AVG_CHI_LOGP),'source','KATE CHI logP','cmpd',pchem.CmpdNumber,'lnb_ref',pchem.ExperimentLNBRef,'prep_lnb_ref',pchem.PrepLNBRef,'date_experiment',pchem.ExperimentDate,'date_record',pchem.ExperimentDate,'date_updated',pchem.DateRecordUpdated); end
                    if ~isempty(pchem.AVG_CHROM_LOGD_PH_7_4),   obj = obj.addPropertyValue('logP','value',str2double(pchem.AVG_CHROM_LOGD_PH_7_4),'source','KATE Chrom LogD pH 7.4','cmpd',pchem.CmpdNumber,'lnb_ref',pchem.ExperimentLNBRef,'prep_lnb_ref',pchem.PrepLNBRef,'date_experiment',pchem.ExperimentDate,'date_record',pchem.ExperimentDate,'date_updated',pchem.DateRecordUpdated); end
                    if ~isempty(pchem.AVG_CHI_LOGD_PH_7_4),     obj = obj.addPropertyValue('logP','value',str2double(pchem.AVG_CHI_LOGD_PH_7_4),'source','KATE CHI LogD pH 7.4','cmpd',pchem.CmpdNumber,'lnb_ref',pchem.ExperimentLNBRef,'prep_lnb_ref',pchem.PrepLNBRef,'date_experiment',pchem.ExperimentDate,'date_record',pchem.ExperimentDate,'date_updated',pchem.DateRecordUpdated); end
                    if ~isempty(pchem.AVG_CHI_IAM),             obj = obj.addPropertyValue('logP','value',str2double(pchem.AVG_CHI_IAM),'source','KATE CHI IAM','cmpd',pchem.CmpdNumber,'lnb_ref',pchem.ExperimentLNBRef,'prep_lnb_ref',pchem.PrepLNBRef,'date_experiment',pchem.ExperimentDate,'date_record',pchem.ExperimentDate,'date_updated',pchem.DateRecordUpdated); end
                end
            end
            
            %CHIN LogP and Acid/Base class
            if ~isempty(CHI_SLST)
                [~,idx] = sort({CHI_SLST.CmpdNumber});
                CHI_SLST = CHI_SLST(idx);
                for i = 1:length(CHI_SLST)
                   chi_s = CHI_SLST(i);
                   if ~isempty([obj.descriptors(strcmpi({obj.descriptors.type},'HbondA')).value])
                       CHIN = max([str2double(chi_s.AVG_CHI_PH_10_5) str2double(chi_s.AVG_CHI_PH_2_0) str2double(chi_s.AVG_CHI_PH_7_4)]);
                       A = obj.descriptors(strcmpi({obj.descriptors.type},'HbondA')).value;
                       CHINlogP = 0.054*CHIN + 1.319*A - 1.877;
                       if ~isnan(CHINlogP)
                           obj = obj.addPropertyValue('logP','value',CHINlogP,'source','KATE CHIN logP','cmpd',chi_s.CmpdNumber,'lnb_ref',chi_s.ExperimentLNBRef,'prep_lnb_ref',chi_s.PrepLNBRef,'date_experiment',chi_s.ExperimentDate,'date_record',chi_s.ExperimentDate,'date_updated',chi_s.DateRecordUpdated);
                       end
                   end
                   obj = obj.addPropertyValue('descriptors','value',chi_s.ACID_OR_BASE,'type','CHI Acid Base class','source','KATE','cmpd',chi_s.CmpdNumber);
                end
            end
            
            %Permeability
            if ~isempty(InVitroPermLST)
                [~,idx] = sort({InVitroPermLST.CmpdNumber});
                InVitroPermLST = InVitroPermLST(idx);
                for i = 1:length(InVitroPermLST)
                    perm = InVitroPermLST(i);
                    if ~isempty(perm.CMPD_PAPP_NM_SEC),         obj = obj.addPropertyValue('perm','Modifier',perm.CMPD_PAPP_MOD,'value',str2double(perm.CMPD_PAPP_NM_SEC),'efflux_ratio',str2double(perm.EFFLUX_RATIO),'pct_recovery',str2double(perm.PCT_RECOVERY),'unit','nm/s','type','Papp','category',sprintf('%s %s pH %s/%s',perm.CellLine,perm.MEAS_DIRECTION,perm.DONOR_PH,perm.RECEIVER_PH),'source','KATE','cmpd',perm.CmpdNumber,'comments',perm.Comments,'study',perm.StudyNumber,'lnb_ref',perm.ExperimentLNBRef,'prep_lnb_ref',perm.PrepLNBRef,'protocol',perm.ProtocolID,'date_experiment',perm.ExperimentDate,'date_record',perm.DateRecordCreated,'date_updated',perm.DateRecordUpdated); end
                    if ~isempty(perm.CMPD_PAPP_INH_NM_SEC),     obj = obj.addPropertyValue('perm','Modifier',perm.CMPD_PAPP_INH_MOD,'value',str2double(perm.CMPD_PAPP_INH_NM_SEC),'efflux_ratio',str2double(perm.EFFLUX_RATIO_INH),'pct_recovery',str2double(perm.PCT_RECOVERY_INH),'unit','nm/s','type','Papp','category',sprintf('%s %s pH %s/%s Inh %s',perm.CellLine,perm.MEAS_DIRECTION,perm.DONOR_PH,perm.RECEIVER_PH,perm.INHIBITOR),'source','KATE','cmpd',perm.CmpdNumber,'comments',perm.Comments,'study',perm.StudyNumber,'lnb_ref',perm.ExperimentLNBRef,'prep_lnb_ref',perm.PrepLNBRef,'protocol',perm.ProtocolID,'date_experiment',perm.ExperimentDate,'date_record',perm.DateRecordCreated,'date_updated',perm.DateRecordUpdated); end
                    if ~isempty(perm.CMPD_PEXACT_NM_SEC),       obj = obj.addPropertyValue('perm','Modifier',perm.CMPD_PEXACT_MOD,'value',str2double(perm.CMPD_PEXACT_NM_SEC),'efflux_ratio',str2double(perm.EFFLUX_RATIO),'pct_recovery',str2double(perm.PCT_RECOVERY),'unit','nm/s','type','Pexact','category',sprintf('%s %s pH %s/%s',perm.CellLine,perm.MEAS_DIRECTION,perm.DONOR_PH,perm.RECEIVER_PH),'source','KATE','cmpd',perm.CmpdNumber,'comments',perm.Comments,'study',perm.StudyNumber,'lnb_ref',perm.ExperimentLNBRef,'prep_lnb_ref',perm.PrepLNBRef,'protocol',perm.ProtocolID,'date_experiment',perm.ExperimentDate,'date_record',perm.DateRecordCreated,'date_updated',perm.DateRecordUpdated); end
                end
            end
            
            %pKa
            if ~isempty(pka)
                for i = 1:length(pka)
                    obj = obj.addPropertyValue('pKa','acidic',pka(i).acid,'basic',pka(i).base,'source','KATE Spectroscopic pKa','study',pka(i).study,'lnb_ref',pka(i).lnb_ref,'prep_lnb_ref',pka(i).prep_lnb_ref,'protocol',pka(i).protocol,'date_experiment',pka(i).date_experiment,'date_record',pka(i).date_record,'date_updated',pka(i).date_updated,'comments',pka(i).comments);
                end
            end
            
            
            %solubility
            if ~isempty(sol)
                sol_clnd = sol(strcmpi({sol.sol_type},'CLND'));
                sol_cad = sol(strcmpi({sol.sol_type},'CAD'));
                sol_fassif = sol(contains({sol.sol_type},'FaSSIF','IgnoreCase',true));
                if ~isempty(sol_clnd)
                    for i = 1:length(sol_clnd)
                        obj = obj.addPropertyValue('sol','cmpd',sol_clnd(i).cmpd,'matrix','Aqueous','pH',sol_clnd(i).pH_measured,'Modifier',sol_clnd(i).sol_mod,'value',sol_clnd(i).sol_mgmL,'unit','mg/mL','source','KATE CLND','prep_lnb_ref',sol_clnd(i).prep_lnb_ref,'lnb_ref',sol_clnd(i).lnb_ref,'protocol',sol_clnd(i).protocol,'date_experiment',sol_clnd(i).date_experiment,'date_record',sol_clnd(i).date_record,'date_updated',sol_clnd(i).date_updated);
                    end
%                     cmpds = unique({sol_clnd.cmpd});
%                     for i = 1:length(cmpds)
%                         sol_iter = sol_clnd(strcmpi({sol_clnd.cmpd},cmpds{i}));
%                         obj = obj.addPropertyValue('sol','cmpd',cmpds{i},'matrix','Aqueous','pH',sol_iter(1).pH_measured,'value',nanmean([sol_iter.sol_mgmL]),'unit','mg/mL','sd',nanstd([sol_iter.sol_mgmL]),'source','KATE CLND');
%                     end
                end
                if ~isempty(sol_cad)
                    for i = 1:length(sol_cad)
                        obj = obj.addPropertyValue('sol','cmpd',sol_cad(i).cmpd,'matrix','Aqueous','pH',sol_cad(i).pH_measured,'Modifier',sol_cad(i).sol_mod,'value',sol_cad(i).sol_mgmL,'unit','mg/mL','source','KATE CAD','prep_lnb_ref',sol_cad(i).prep_lnb_ref,'lnb_ref',sol_cad(i).lnb_ref,'protocol',sol_cad(i).protocol,'date_experiment',sol_cad(i).date_experiment,'date_record',sol_cad(i).date_record,'date_updated',sol_cad(i).date_updated);
                    end
%                     cmpds = unique({sol_cad.cmpd});
%                     for i = 1:length(cmpds)
%                         sol_iter = sol_cad(strcmpi({sol_cad.cmpd},cmpds{i}));
%                         obj = obj.addPropertyValue('sol','cmpd',cmpds{i},'matrix','Aqueous','pH',sol_iter(1).pH_measured,'value',nanmean([sol_iter.sol_mgmL]),'unit','mg/mL','sd',nanstd([sol_iter.sol_mgmL]),'source','KATE CAD');
%                     end
                end
                if ~isempty(sol_fassif)
                    for i = 1:length(sol_fassif)
                        obj = obj.addPropertyValue('sol','cmpd',sol_fassif(i).cmpd,'matrix',sol_fassif(i).solvent,'pH',sol_fassif(i).pH_measured,'Modifier',sol_fassif(i).sol_mod,'value',sol_fassif(i).sol_mgmL,'unit','mg/mL','source','KATE FaSSIF HT','prep_lnb_ref',sol_fassif(i).prep_lnb_ref,'lnb_ref',sol_fassif(i).lnb_ref,'protocol',sol_fassif(i).protocol,'date_experiment',sol_fassif(i).date_experiment,'date_record',sol_fassif(i).date_record,'date_updated',sol_fassif(i).date_updated);
                    end
%                     cmpds = unique({sol_fassif.cmpd});
%                     for i = 1:length(cmpds)
%                         sol_iter = sol_fassif(strcmpi({sol_fassif.cmpd},cmpds{i}));
%                         obj = obj.addPropertyValue('sol','cmpd',cmpds{i},'matrix',sol_iter(1).solvent,'pH',sol_iter(1).pH_measured,'value',nanmean([sol_iter.sol_mgmL]),'unit','mg/mL','sd',nanstd([sol_iter.sol_mgmL]),'source','KATE FaSSIF HT');
%                     end
                end
            end
            
            %Protein Binding
            if ~isempty(ProteinBindingLST)
                bnd = struct('cmpd',    {ProteinBindingLST.CmpdNumber},...
                             'species', [],...
                             'strain',  [],...
                             'sex',     [],...
                             'matrix',  {ProteinBindingLST.MATRIX},...
                             'method',  {ProteinBindingLST.METHOD},...
                             'fu_mod',  {ProteinBindingLST.CMPD_PCT_BINDING_MOD},...
                             'fu_pct',  num2cell(str2double({ProteinBindingLST.FU_PERCENT})./100),...
                             'conc',    num2cell(str2double({ProteinBindingLST.CMPD_CONC})),...
                             'conc_units',{ProteinBindingLST.CMPD_CONC_UNITS},...
                             'n',       num2cell(str2double({ProteinBindingLST.REPLICATES})),...
                             'std',     num2cell(str2double({ProteinBindingLST.FU_SD})),...
                             'study',   {ProteinBindingLST.StudyNumber},...
                             'date_experiment', {ProteinBindingLST.ExperimentDate},...
                             'date_record', {ProteinBindingLST.DateRecordCreated},...
                             'date_updated',{ProteinBindingLST.DateRecordUpdated},...
                             'lnb_ref',     {ProteinBindingLST.ExperimentLNBRef},...
                             'prep_lnb_ref',{ProteinBindingLST.PrepLNBRef},...
                             'protocol',    {ProteinBindingLST.ProtocolID},...
                             'comments',{ProteinBindingLST.Comments});
                for i = 1:length(ProteinBindingLST)
                    bnd(i).species = lower(ProteinBindingLST(i).PBanimal.Species);
                    bnd(i).strain  = lower(ProteinBindingLST(i).PBanimal.Strain);
                    bnd(i).sex     = lower(ProteinBindingLST(i).PBanimal.Gender);
                end
                
                sp = unique(lower({bnd.species}));
                ma = unique(lower({bnd.matrix}));
                cn = unique({bnd.cmpd});
                for i = 1:length(sp)
                    for j = 1:length(ma)
                        for k = 1:length(cn)
                            details = bnd(strcmpi({bnd.species},sp{i}) & strcmpi({bnd.matrix},ma{j}) & strcmpi({bnd.cmpd},cn{k}));
                            if ~isempty(details)
                                for l = 1:length(details)
                                    det(1).strain = details(l).strain;
                                    det(1).sex = details(l).sex;
                                    det(1).fu_mod = details(l).fu_mod;
                                    if isnan(details(l).conc)
                                        det(1).conc = [];
                                    else
                                        det(1).conc = details(l).conc;
                                    end
                                    det(1).conc_units = details(l).conc_units;
                                    if isnan(details(l).n)
                                        det(1).n = [];
                                    else
                                        det(1).n = details(l).n;
                                    end
                                    if isnan(details(l).std)
                                        det(1).std = [];
                                    else
                                        det(1).std = details(l).std;
                                    end
                                    obj = obj.addPropertyValue('bind','cmpd',cn{k},'matrix',ma{j},'fu',details(l).fu_pct,'species',sp{i},'details',det,'source','KATE','method',details(l).method,'comments',details(l).comments,'date_experiment',details(l).date_experiment,'date_record',details(l).date_record,'date_updated',details(l).date_updated,'lnb_ref',details(l).lnb_ref,'prep_lnb_ref',details(l).prep_lnb_ref,'protocol',details(l).protocol);
                                end
                            end
                        end
                    end
                end
            end
            
            %fup
            bnd = obj.bind(strcmpi({obj.bind.source},'kate'));
            sp = unique({bnd.species});
            cn = unique({bnd.cmpd});
            for i = 1:length(sp)
                for j = 1:length(cn)
                    idx = find(strcmp({bnd.species},sp{i}) & strcmp({bnd.cmpd},cn{j}) & strcmpi({bnd.matrix},'plasma'));
                    if ~isempty(idx)
                        for k = 1:length(bnd(idx))
                            obj = obj.addPropertyValue('fup','cmpd',cn{j},'value',bnd(idx(k)).fu,'species',sp{i},'details',bnd(idx(k)).details,'source','KATE','method',bnd(idx(k)).method,'comments',bnd(idx(k)).comments,'study',bnd(idx(k)).study,'lnb_ref',bnd(idx(k)).lnb_ref,'prep_lnb_ref',bnd(idx(k)).prep_lnb_ref,'protocol',bnd(idx(k)).protocol,'date_experiment',bnd(idx(k)).date_experiment,'date_record',bnd(idx(k)).date_record,'date_updated',bnd(idx(k)).date_updated);
                        end
                    end
                end
            end
            
            % HSA, AGP, and Max Drug Efficiency
            if ~isempty(PhysChemSumLST)
                [~,idx] = sort({PhysChemSumLST.CmpdNumber});
                PhysChemSumLST = PhysChemSumLST(idx);
                for i = 1:length(PhysChemSumLST)
                    pchem = PhysChemSumLST(i);
                    if ~isempty(pchem.AGP_PCT_BIND_MEAN)
                        obj = obj.addPropertyValue('bind','cmpd',pchem.CmpdNumber,'matrix','AGP','fu',(100 - (str2double(pchem.AGP_PCT_BIND_MEAN)))/100,'source','KATE','method','HT AGP','species','human');
                    end
                    if ~isempty(pchem.HSA_PCT_BIND_MEAN)
                        obj = obj.addPropertyValue('bind','cmpd',pchem.CmpdNumber,'matrix','HSA','fu',(100 - (str2double(pchem.HSA_PCT_BIND_MEAN)))/100,'source','KATE','method','HT HSA','species','human');
                    end
                    if ~isempty(pchem.HPLC_DRUG_EFF_MAX_PCT)
                        obj = obj.addPropertyValue('logP','cmpd',pchem.CmpdNumber,'value',str2double(pchem.HPLC_DRUG_EFF_MAX_PCT),'source','KATE_HPLC_DRUG_MAX_EFF_PCT');
                    end
                end
            end
            
            %B:P ratio
            if ~isempty(BloodPlasmaPartitionLST)
                b_p = BloodPlasmaPartitionLST;
                for i = 1:length(b_p)
                    details = struct('Species',[],'Strain',[],'Sex',[],'Concentration',[],'Units',[],'BP_Ratio',[]);
                    details(1).Concentration = str2double(b_p(i).CONC1_NG_ML);   details(1).BP_Ratio = str2double(b_p(i).B_P_RATIO1);
                    details(2).Concentration = str2double(b_p(i).CONC2_NG_ML);   details(2).BP_Ratio = str2double(b_p(i).B_P_RATIO2);
                    details(3).Concentration = str2double(b_p(i).CONC3_NG_ML);   details(3).BP_Ratio = str2double(b_p(i).B_P_RATIO3);
                    [details.species] = deal(lower([b_p(i).BPPAnimal.Species]));
                    [details.strain]  = deal(lower([b_p(i).BPPAnimal.Strain]));
                    [details.sex]     = deal(b_p(i).BPPAnimal.Gender);
                    if std([details.BP_Ratio]) > (0.25 * mean([details.BP_Ratio])), flag = 'Standard deviation of B:P Ratio greater than 25% of the mean. This may indicate concentration dependence'; else, flag = ''; end
                    obj = obj.addPropertyValue('bp','cmpd',b_p(i).CmpdNumber,'value',nanmean([details.BP_Ratio]),'std',nanstd([details.BP_Ratio]),'n',str2double(b_p(i).BPPAnimal.NumAnimals),'species',lower(b_p(i).BPPAnimal.Species),'details',details,'source','KATE','flag',flag,'comments',b_p(i).Comments,'study',b_p(i).StudyNumber,'lnb_ref',b_p(i).ExperimentLNBRef,'prep_lnb_ref',b_p(i).PrepLNBRef,'protocol',b_p(i).ProtocolID,'date_experiment',b_p(i).ExperimentDate,'date_record',b_p(i).DateRecordCreated,'date_updated',b_p(i).DateRecordUpdated);
                end
            end
            
            %Intrinsic clearance
            if ~isempty(InVitroMetStabLST)
                ms = InVitroMetStabLST;
                for i = 1:length(ms)
                    if ~isnan(str2double(ms(i).CL_ML_MIN_G_TISSUE))
                        details = struct('species',lower(ms(i).MetStabilityAnimal.Species),'strain',lower(ms(i).MetStabilityAnimal.Strain),'sex',lower(ms(i).MetStabilityAnimal.Gender),'CmpdConc',ms(i).COMD_CONC_UM,'CmpdConcUnit','uM','SystemConc',ms(i).INVITRO_SYSTEM_CONC,'SystemUnit',ms(i).INVITRO_SYSTEM_UNIT);
                        obj = obj.addPropertyValue('met','cmpd',ms(i).CmpdNumber,'Modifier',ms(i).CL_ML_MIN_G_TISSUE_MOD,'CL',str2double(ms(i).CL_ML_MIN_G_TISSUE),'unit','mL/min/g','species',lower(ms(i).MetStabilityAnimal.Species),'system',lower(ms(i).IN_VITRO_SYSTEM),'details',details,'source','KATE','study',ms(i).StudyNumber,'lnb_ref',ms(i).ExperimentLNBRef,'prep_lnb_ref',ms(i).PrepLNBRef,'protocol',ms(i).ProtocolID,'date_experiment',ms(i).ExperimentDate,'date_record',ms(i).DateRecordCreated,'date_updated',ms(i).DateRecordUpdated,'comments',ms(i).Comments);
                    end
                end
            end
        end
        
        function obj = getKATEPKProperties(obj)
            url 	= 'https://ivivtapps-prod.gsk.com/KATESearch/KATESearch.svc/REST/';
            format  = 'json';
            opt     = weboptions('Timeout',30,'CertificateFilename','');
            
            timer = tic;
                InVivoPharmacokineticsLST = webread([url format '/' 'GetInVivoPharmacokineticsResults' '?key=0&cmpd=' obj.name],opt);
                if isempty(InVivoPharmacokineticsLST), InVivoPharmacokineticsLST = webread([url format '/' 'GetInVivoPharmacokineticsResults' '?key=0&cmpd=' obj.gcn],opt); end
            timer = toc(timer);
            if isempty(InVivoPharmacokineticsLST)
                obj = obj.addPropertyValue('fxn_status','service','getKATEPKProperties','time_s',timer,'ReturnedValues',0,'comments',sprintf('No PK data found for %s',obj.name));
            else
                obj = obj.addPropertyValue('fxn_status','service','getKATEPKProperties','time_s',timer,'ReturnedValues',1);
            end
            
            if ~isempty(InVivoPharmacokineticsLST)
                pk = InVivoPharmacokineticsLST;
                for i = 1:length(pk)
                    obj = obj.addPropertyValue('pk_kate','cmpd',            pk(i).CmpdNumber,...
                                                         'study',           pk(i).StudyNumber,...
                                                         'lnb_ref',         pk(i).ExperimentLNBRef,...
                                                         'prep_lnb_ref',    pk(i).PrepLNBRef,...
                                                         'study_design',    lower(pk(i).STUDY_DESIGN),...
                                                         'species',         lower(pk(i).PKanimal.Species),...
                                                         'strain',          lower(pk(i).PKanimal.Strain),...
                                                         'sex',             lower(pk(i).PKanimal.Gender),...
                                                         'animal_id',       pk(i).PKanimal.AnimalID,...
                                                         'n',               str2double(pk(i).TOT_NBR_ANIMALS),...
                                                         'dose_mgkg',       str2double(pk(i).DOSE_MG_KG),...
                                                         'route',           pk(i).ROUTE_ADMIN,...
                                                         'formulation',     pk(i).FORMULATION,...
                                                         'matrix',          lower(pk(i).SAMPLE_MATRIX),...
                                                         'F_mod',           pk(i).F_PCT_MOD,...
                                                         'F',               str2double(pk(i).F_PCT),...
                                                         'F_std',           str2double(pk(i).F_PCT_SD),...
                                                         'Cmax_mod',        pk(i).CMAX_MOD,...
                                                         'Cmax_ngml',       str2double(pk(i).CMAX_NG_ML),...
                                                         'Cmax_std',        str2double(pk(i).CMAX_SD),...
                                                         'Tmax_mod',        pk(i).TMAX_MOD,...
                                                         'Tmax_hr',         str2double(pk(i).TMAX_HR),...
                                                         'Tmax_std',        str2double(pk(i).TMAX_SD),...
                                                         'AUC_Inf_mod',     pk(i).AUC_0INF_MOD,...
                                                         'AUC_Inf_nghrmL',  str2double(pk(i).AUC_0INF_NG_HR_ML),...
                                                         'AUC_Inf_std',     str2double(pk(i).AUC_0INF_SD),...
                                                         'AUC_t_mod',       pk(i).AUC_0T_MOD,...
                                                         'AUC_t_nghrmL',    str2double(pk(i).AUC_0T_NG_HR_ML),...
                                                         'AUC_t_std',       str2double(pk(i).AUC_0T_SD),...
                                                         'Clast_mod',       pk(i).C_LAST_MOD,...
                                                         'Clast_ngml',      str2double(pk(i).C_LAST_NG_ML),...
                                                         'Clast_std',       str2double(pk(i).C_LAST_SD),...
                                                         'Tlast_mod',       pk(i).T_LAST_MOD,...
                                                         'Tlast_hr',        str2double(pk(i).T_LAST_HR),...
                                                         'TBdy_CL_mod',     pk(i).TOT_BDY_CLR_MOD,...
                                                         'TBdy_CL_mlmin_kg',str2double(pk(i).TOT_BDY_CLR_ML_MIN_KG),...
                                                         'TBdy_CL_std',     str2double(pk(i).TOT_BDY_CLR_SD),...
                                                         'Vdss_mod',        pk(i).VDSS_MOD,...
                                                         'Vdss_Lkg',        str2double(pk(i).VDSS_L_KG),...
                                                         'Vdss_std',        str2double(pk(i).VDSS_SD),...
                                                         'T1_2_mod',        pk(i).T1_2_MOD,...
                                                         'T1_2_hr',         str2double(pk(i).T1_2_HR),...
                                                         'T1_2_std',        str2double(pk(i).T1_2_SD),...
                                                         'MRT_mod',         pk(i).MRT_MOD,...
                                                         'MRT_hr',          str2double(pk(i).MRT_HR),...
                                                         'MRT_std',         str2double(pk(i).MRT_SD),...
                                                         'FA_HPV_mod',      pk(i).FA_HPV_MOD,...
                                                         'FA_HPV',          str2double(pk(i).FA_HPV),...
                                                         'FA_HPV_std',      str2double(pk(i).FA_HPV_SD),...
                                                         'project',         pk(i).ProjectName,...
                                                         'protocol_ID',     pk(i).ProtocolID,...
                                                         'pk_start_date',   pk(i).ExperimentDate,...
                                                         'pk_date_record',  pk(i).DateRecordCreated,...
                                                         'pk_date_updated', pk(i).DateRecordUpdated,...
                                                         'comments',        pk(i).Comments,...
                                                         'source',          'KATE');
                    flds = fields(obj.pk_kate(end));                
                    for j = 1:length(flds)
                        if prod(isnan(obj.pk_kate(i).(flds{j}))), obj.pk_kate(i).(flds{j}) = []; end
                    end
                end
                [~,idx] = sort([obj.pk_kate.dose_mgkg]);obj.pk_kate = obj.pk_kate(idx);
                [~,idx] = sort({obj.pk_kate.species});  obj.pk_kate = obj.pk_kate(idx);
                [~,idx] = sort({obj.pk_kate.cmpd});     obj.pk_kate = obj.pk_kate(idx);
            end 
            
        end
        
        function obj = getDatamartProperties(obj)
        end
    end
end

