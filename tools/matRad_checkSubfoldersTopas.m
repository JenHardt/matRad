function out = matRad_checkSubfoldersTopas(folder)

dirinfo = dir(folder);
counter = 0;
for d = 1:length(dirinfo)
    if length(dirinfo(d).name) > 3
        counter = counter + 1;
        out(counter) = length(dir([dirinfo(d).folder filesep dirinfo(d).name filesep '*.bin'])) > 3;
    end
end

end