function s = getInfoSingleVirtualCompound(vcmpds, idx)
s = struct();
for ii = string(fieldnames(vcmpds)')
    temp = vcmpds.(ii);
    if length(temp) > 1
        s.(ii) = temp(idx);
    else
        s.(ii) = temp(1);
    end
end
end