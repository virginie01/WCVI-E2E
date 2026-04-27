function BioIn = loadBiologicalForcing(Forcingdir, LBdir)
%LOADBIOLOGICALFORCING Load biological forcing datasets from .mat files

files = {
    'o2_input'
    'fisheriesTS_input'
};

for i = 1:numel(files)
    fname = files{i};
    S = load(fullfile(Forcingdir, [fname '.mat']));
    BioIn.(fname) = S.(fname);
end


files = dir(fullfile(LBdir, '*input.mat'));

for i = 1:length(files)
    tmp = load(files(i).name);
    fld = fieldnames(tmp);
    BioIn.LB(i).data = tmp.(fld{1});
    BioIn.LB(i).name = fld{1};
end

end
