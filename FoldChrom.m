function[Mod,uRt] = FoldChrom(timeVector,modTime,trace)

%%
Options.modTime =modTime;
wU = Options.modTime;
Mod = zeros(length(timeVector),1);
Mod(1) = 1;
for n = 2:length(timeVector)
if timeVector(n) > wU & trace(n) == 1
Mod(n) = 1;
wU = wU + Options.modTime;
% wU + Options.modTime;
end 
end 
Mod = cumsum(Mod);
medianTime = median(diff(unique(timeVector)));

[~,~,uRt] =unique(timeVector);% round([medianTime;diff(timeVector)]/medianTime);
% old version
% uRt = round([medianTime;diff(timeVector)]/medianTime);
% uRt = cumsum(uRt);
% uRt = round_odd(uRt);

[~,a,ind] = unique(Mod); % I'm here
% a = [1;a(1:end-1)];
uRt =uRt-(uRt(a(ind)))+1;
% uRt(ind == 1) = uRt(ind == 1)+1;

end 
function S = round_odd(S)
% round to nearest odd integer.
idx = mod(S,2)<1 & S>0;
S = floor(S);
S(idx) = S(idx)+1;
end 