function out = getSelectDPSmodelOutputs(smi,model,varargin)
    %getSelectDPSmodelOutputs returns the requested outputs of any model
    %loaded on the DPS (Derived Property Server)
    %   Inputs
    %       1) smi - SMILES string or gcn of the requested compound. If
    %       using gcn, change the searchType to gcn (see name-value
    %       arguments below)
    %       2) model - Model from which outputs are requested. This can be
    %       strictly a model name (e.g., CLiv_rat_v1) which returns all
    %       properties, or a model name with a property name separated by a
    %       . (e.g., CLiv_rat_v1.CL_class) which only returns the value for
    %       property CL_class. This can be a char or a cell array of
    %       multiple models/model.Property
    %      
    %       3) Name-Value pairs
    %           a) searchType - specifies the search input. Valid options
    %           are smiles or gcn, default is smiles. If searching by gcn,
    %           be sure to include "searchType,gcn" in the function call
    %
    %   Output
    %       1) Structure of all requested models and properties
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    searchType = 'smiles';
    for i = 1:length(varargin)
        switch lower(varargin{i})
            case 'searchtype', searchType = varargin{i+1};
        end
    end
    switch lower(searchType)
        case 'gcn'
            gcn = smi;
            smi = smtb.webserv.GetSMILESFromCLrest(gcn);
            if isempty(smi)
                error('SMILES string not found for compound %s',gcn)
            end
        case 'smiles' % do nothing
        otherwise, error('Invalid search type %s. Please choose from smiles or gcn',searchType)
    end
    
    if iscell(model)
        for i = 1:length(model) - 1
            model{i} = [model{i} ';'];
        end
        model = [model{:}];
    end
    
    smi_str = '';
    if iscell(smi)
        for i = 1:length(smi)
            if i == length(smi)
                smi_str = cat(2,smi_str,smi{i});
            else
                smi_str = cat(2,smi_str,smi{i},';');
            end
        end
    else
        smi_str = smi;
    end
    
    smi_str = strrep(smi_str,'#','%23');
    
    options = weboptions;
    options.CertificateFilename = '';
    options.Timeout = 180;
    
%     url = sprintf('https://ivivtapps-prod.gsk.com/SMTBToolboxService/SMTBToolboxService.svc/REST/json/GetSelectCalculatedDerivedPropertiesDPS?props=%s&smi=%s',model,smi);
    url = sprintf('https://ivivtapps-prod.gsk.com/SMTBToolboxService/SMTBToolboxService.svc/REST/json/GetSelectCalculatedDerivedPropertiesDPS?props=%s&smi=%s',model,smi_str);
    
    out = webread(url,options);
    
%     out_read = webread(url,options);
%     out = out_read;
%     
%     out = struct.empty;
%     
%     if ~isempty(out_read)
%         props = {out_read.PropertyName};
% 
%         if contains(model,'.')
%             spl = strsplit(model,'.');
%             model = spl{1};
%         end
%         count = 0;
%         for i = 1:length(props)
%             if contains(props{i},model)
%                 spl = strsplit(props{i},'.');
%                 propi = spl{2};
%             else
%                 propi = props{i};
%             end
%             count = count+1;
%             out(1).(propi) = out_read(strcmpi({out_read.PropertyName},props{i})).PropertyValue;
%         end
%     end
    
end

