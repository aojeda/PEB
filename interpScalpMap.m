function [M, Sout] = interpScalpMap(chanlocsIn, chanlocsOut, Sin)
% Computes the linear mapping M: In chanlocs space -> Out chanlocs space,
% where chanlocsIn and chanlocsOut are EEGLAB channel location structures
% that correspond to the In and out spaces respectively.   

xyzOut = [[chanlocsOut.X]' [chanlocsOut.Y]' [chanlocsOut.Z]'];
xyzIn =  [[chanlocsIn.X]'  [chanlocsIn.Y]'  [chanlocsIn.Z]'];
labelsIn  = {chanlocsIn.labels};
labelsOut = {chanlocsOut.labels};
Nout = size(xyzOut,1);
Nin = size(xyzIn,1);
M = zeros(Nout,Nin);
Nneig = 4;
for k=1:Nout
    I = ismember(labelsIn,labelsOut{k});
    if any(I)
        M(k,I) = 1;
    else
        d = sum(bsxfun(@minus,xyzIn,xyzOut(k,:)).^2,2);
        [~,sorting] = sort(d);
        m = 1./d(sorting(1:Nneig));
        m = m/sum(m);
        M(k,sorting(1:Nneig)) = m;
    end
end
if nargin > 2
    Sout = M*Sin;
else
    Sout = [];
end
end