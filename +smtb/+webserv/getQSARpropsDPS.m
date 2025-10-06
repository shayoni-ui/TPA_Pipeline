function profile = getQSARpropsDPS(input,varargin)
    %getQSARpropsDPS
    %   Detailed explanation goes here
    
%     url = 'https://ivivtapps-test.gsk.com/SMTBToolboxService/SMTBToolboxService.svc/REST/json/';
    url = 'https://ivivtapps-prod.gsk.com/SMTBToolboxService/SMTBToolboxService.svc/REST/json/';
    db_service  = 'GetQSARPropertiesByCmpd';
%     db_service  = 'addservice'; add input option to run models as opposed
%     to first searching the cache.
    
    options = weboptions;
    options.Timeout = 180;
    options.MediaType = 'application/x-www-form-urlencoded';
    options.CertificateFilename = '';
    
    if ischar(input), input = {input}; end
    if ~isempty(varargin)
        if ischar(varargin{1}), varargin{1} = varargin(1); end
    end
    if size(input,2)>size(input,1), input = input'; end

    gcn = cell(size(input,1),1);
    smiles = cell(size(input,1),1);

    for i = 1:length(input)
       if smtb.webserv.isGSKcn(input{i})==1
           gcn{i} = input{i};
           smiles{i} = '';  % ADMET service will add the SMILES from CODS
       else
           smiles{i} = input{i};
           switch length(varargin)
               case 0, gcn{i} = ''; % do nothing
               case 1
                   if smtb.webserv.isGSKcn(varargin{1}{i})==1
                       gcn{i} = varargin{1}{i};
                       if any(contains(smiles{i},'.'))||any(contains(smiles{i},'error','IgnoreCase',true))
                           smiles{i} = '';
                           gcn{i} = varargin{1}{i};
                           warning('SMILES string is for salt form.  Will attempt to retrieve parent SMILES.')
                       end
                   else
                      gcn{i} = '';
                      if any(contains(smiles{i},'.'))||any(contains(smiles{i},'error','IgnoreCase',true))
                           smiles{i} = strrep(smiles{i},'.','');
                           gcn{i} = varargin{1}{i};
                           warning('SMILES string is for salt form.  Will attempt to retrieve parent SMILES.')
                       end
                   end
               otherwise
                   warning('Optional input is reserved for GSK compound number.  Additional input is being ignored')
                   gcn{i} = '';
           end
           if any(contains(smiles{i},'.'))||any(contains(smiles{i},'error','IgnoreCase',true))
               smiles{i} = strrep(smiles{i},'.','');
               warning('SMILES string is for salt form.  Calculated properties will be for salt form.  Provide GSK compound number as second argument to retrieve parent SMILES')
           end
       end
    end
    
    if size(gcn,2)>size(gcn,1), gcn = gcn'; end
    if size(smiles,2)>size(smiles,1), smiles = smiles'; end
    
    profile = cell(length(gcn),1);
    for i = 1:length(gcn)
        run_url = [url db_service];
        p = webread(sprintf('%s?cmpd=%s&smi=%s',run_url,gcn{i},smiles{i}),options);
        profile{i} = p;
    end
        
    profile = [profile{:}];
    
end

