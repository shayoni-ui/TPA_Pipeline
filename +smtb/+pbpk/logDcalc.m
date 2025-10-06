function logD = logDcalc( varargin )
%LOGDCALC Calculates logD from logP, pH, and pKa's using the GastroPlus
%empirical logD correction factors
%
%   logD = logDcalc(logP,pH,Acidic pKas,Basic pKas)
%   logD = logDcalc(drugprops,pH)

if isa(varargin{1},'DrugProps')
    logP  = getPreferredValue(varargin{1},'logP');
    pKa_a = getPreferredValue(varargin{1},'pKa_a');
    pKa_b = getPreferredValue(varargin{1},'pKa_b');
    pH    = smtb.useful.typecheck(varargin{2},'double','pH');
    logD  = logDcalc(logP,pH,pKa_a,pKa_b);
    return
else
    logP  = smtb.useful.typecheck(varargin{1},'double','logP');
    pH    = smtb.useful.typecheck(varargin{2},'double','pH');
    pKa_a = smtb.useful.typecheck(varargin{3},'double','Acidic pKa');
    pKa_b = smtb.useful.typecheck(varargin{4},'double','Basic pKa');
end

PaPn = 3.36; %logP(neutral)-logP(anion)
PcPn = 3.06; %logP(neutral)-logP(cation)
PzPn = 2.44; %logP(neutral)-logP(zwitterion)

f     = smtb.pbpk.ionfrac(pKa_a,pKa_b,pH);

logD  = log10(f.neutral      * 10^logP +...
              f.cationic     * 10^(logP-PcPn) + ...
              f.anionic      * 10^(logP-PaPn) + ...
              f.zwitterionic * 10^(logP-PzPn));

end