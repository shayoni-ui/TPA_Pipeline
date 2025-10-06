function profile = getADMETProps(input,varargin)
%GETADMETPROPS
%   Detailed explanation goes here

url     = 'https://ivivtapps-prod.gsk.com/ADMETPropertiesService/ADMETPropertiesService.svc/REST/json'; % Compatible with JSON or XML, but designed for JSON use
% url     = 'https://ivivtapps-test.gsk.com/ADMETPropertiesService/ADMETPropertiesService.svc/REST/json'; % Compatible with JSON or XML, but designed for JSON use

% service = 'BatchProcessADMETCompoundDynamicByStructure';
service = 'BuildMasterADMETCompoundProfile';
url = [url '/' service];

options = weboptions;
options.Timeout = 60;
options.MediaType = 'application/x-www-form-urlencoded';
options.CertificateFilename = '';

if ischar(input), input = {input}; end
if ~isempty(varargin)
    if length(varargin) == 1
        if ischar(varargin{1})
            varargin{1} = varargin(1);
        elseif isa(varargin{1},'double')
            options.Timeout = min(300,varargin{1});
            varargin = {};
        end
    elseif length(varargin) == 2
        options.Timeout = min(300,varargin{2});
        if ischar(varargin{1})
            varargin = {varargin(1)};
        else
            varargin = varargin(1);
        end
    end
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

if length(input) <= 100
    Compound = gcn;
    Structure = smiles;
    StructureFormat = cell(size(Structure,1),1); [StructureFormat{:}] = deal('SMI');
    jsonInput = jsonencode(table(Compound, Structure, StructureFormat));
    profile = webwrite(url,jsonInput,options);
else
    for i = 1:ceil(length(input)/100)
        if i==ceil(length(input)/100)
            Compound = gcn(((i*100)-99):length(input));
            Structure = smiles(((i*100)-99):length(input));
            StructureFormat = cell(size(Structure,1),1); [StructureFormat{:}] = deal('SMI');
            jsonInput = jsonencode(table(Compound, Structure, StructureFormat));
            profile(((i*100)-99):length(input)) = webwrite(url,jsonInput,options);
        else
            Compound = gcn(((i*100)-99):i*100);
            Structure = smiles(((i*100)-99):i*100);
            StructureFormat = cell(size(Structure,1),1); [StructureFormat{:}] = deal('SMI');
            jsonInput = jsonencode(table(Compound, Structure, StructureFormat));
            profile(((i*100)-99):i*100) = webwrite(url,jsonInput,options);
        end
    end
end 
% else
%     service = 'ProcessADMETCompoundDynamicByStructure';
%     strct = char(smiles);
%     cmpd = '';
%     fmt  = 'SMI';
%     
% %     switch length(varargin)
% %         case 0
% %         case 1, format = varargin{1};
% %         otherwise, error('Too many input arguments, please see help getADMETProps');
% %     end
% % 
% %     switch lower(format)
% %         case 'gcn',    cmpd = gcn; strct = '';  fmt = 'SMI';
% %         case 'smiles', cmpd = '';  strct = gcn; fmt = 'SMI';
% %         case 'mol',    cmpd = '';  strct = gcn; fmt = 'MOL';
% %         case 'sdf',    cmpd = '';  strct = gcn; fmt = 'SDF';
% %         otherwise, error('Unrecognized input format %s',format)
% %     end
% 
%     url = [url '/' service '?cmpd=' cmpd '&struct=' strct '&fmt=' fmt];
%     url = strrep(url_0,'#','%23');  %Remove # from url with equivalent %23 character
% 
%     % Turn off certificate as this was causing errors. Considered low risk as
%     % this is an intranet service
%     opt     = weboptions('Timeout',60,'CertificateFilename','');
%     profile = webread(url,opt);
% end

end

