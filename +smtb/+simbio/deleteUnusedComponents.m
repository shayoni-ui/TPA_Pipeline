function deleteUnusedComponents(mm)
% delete unsused components in the simbiology model
%

for jj = {'parameter', 'species', 'compartment'}
    x=sbioselect(mm,'type',jj{1});
    l=zeros(size(x));
    for k=1:length(x)
        l(k)=length(findUsages(x(k)));
    end
    x(l==0).delete;
end

end