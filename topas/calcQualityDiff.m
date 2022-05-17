function qiComp = calcQualityDiff(qi1,qi2)

fnames = fieldnames(qi1);
qiComp = zeros(length(qi1),length(fnames)-1);
for f = 2:length(fnames)-8
    for j = 1:length({qi1.(fnames{f})})
        x = (qi1(j).(fnames{f})-qi2(j).(fnames{f}))./qi1(j).(fnames{f})*100;
        if isscalar(x)
            qiComp(j,f) = x;
        end
    end
end

end

