function sortBasedata(pathToBasedata)

load(pathToBasedata,'machine')

fnames = lower(fieldnames(machine.data));
alphaBeta = contains(fnames,{'alpha','beta','LET','AbsNb'});
fnamesAB = fnames;
fnamesAB(alphaBeta) = deal({'zzzz'});
[~,idx] = sort(fnamesAB);

machine.data = orderfields(machine.data, idx);

save(pathToBasedata,'machine')

end