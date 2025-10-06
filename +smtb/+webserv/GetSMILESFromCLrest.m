function GetSMILESFromCLResult = GetSMILESFromCLrest(sCmpdNumber)
%GetSMILESFromCL(obj,sCmpdNumber)
%
%   Method for retrieving SMILE string from Chemistry Lookup Service based on compound number.
%   
%     Input:
%       sCmpdNumber = (string)
%   
%     Output:
%       GetSMILESFromCLResult = (string)

% Build up the argument lists.
% values = { ...
%    sCmpdNumber, ...
%    };
% names = { ...
%    'sCmpdNumber', ...
%    };
% types = { ...
%    '{http://www.w3.org/2001/XMLSchema}string', ...
%    };
% 
% % Create the message, make the call, and convert the response into a variable.
% soapMessage = createSoapMessage( ...
%     'http://dmpkwikidev:8086/MTBModelBrowser', ...
%     'GetSMILESFromCL', ...
%     values,names,types,'document');
% response = callSoapService( ...
%     obj.endpoint, ...
%     'http://dmpkwikidev:8086/MTBModelBrowser/GetSMILESFromCL', ...
%     soapMessage);
% GetSMILESFromCLResult = parseSoapResponse(response);


options = weboptions;
options.CertificateFilename = '';
options.Timeout = 60;

gcn = smtb.useful.typecheck(sCmpdNumber,'char','gcn');

url_cmpd = sprintf('https://ivivtapps-prod.gsk.com/SMTBToolboxService/SMTBToolboxService.svc/REST/json/GetSMILESFromCL?cmpd=%s',gcn);

GetSMILESFromCLResult = webread(url_cmpd,options);