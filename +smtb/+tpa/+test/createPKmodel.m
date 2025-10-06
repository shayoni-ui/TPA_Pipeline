function [m, mexp] = createPKmodel(nva)
arguments
    nva.exportFlag=false;
end
    newcomp.name = 'TargetTissue';
    newcomp.capacity = 1;
    newcomp.Q = 100;
    [m, mexp] = smtb.tpa.selectPKmodel('PBPK', ...
        'newCompartment',newcomp,'exportFlag',nva.exportFlag);
end