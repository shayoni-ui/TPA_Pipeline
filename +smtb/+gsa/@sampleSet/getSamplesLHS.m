function samples = getSamplesLHS(obj)
    if obj.Seed
        rng(obj.Seed);
    end
    samples = lhsdesign(obj.NS, obj.NP);
    %sampleDist = cell(1, obj.NP);
    %[sampleDist{:}] = deal(makedist('Uniform'));
end