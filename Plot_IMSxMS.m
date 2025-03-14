function [Location,mz] = Plot_IMSxMS(mzroi,MSroi,Mod,uRt,Options,CompoundNumber,trace)
fig = figure;
% fig.Position = [-2167 1232 958 1066];
col = 'rb';
indTrace=unique(uRt(trace == 1));

for T = 1: Options.NumTrace
    subplot(2,1,T)
    dtInd = (Mod == Options.Location(CompoundNumber,1) & uRt == Options.Location(CompoundNumber,2)+T-1);
    StepSize = sum(dtInd);
    
%     dtInd = dtInd(:);
nVec = -Options.SumNumScan:Options.SumNumScan;
Z = zeros(size(MSroi,1),sum(dtInd));
dtInd = find( dtInd);
for n = 1:length(nVec)
     dtInt_tmp = dtInd+nVec(n) * StepSize * 2;
    Z = full(MSroi(:,dtInt_tmp))+Z;
end 
%     Z = accumarray([repmat(1:size(Z,1),1,Options.SumNumScan*2+1),[repmat(1:StepSize,1,Options.SumNumScan*2+1)])
    indMZ = sum(Z>0,2)>0;
    Z = Z(indMZ,:);
    [X,Y] = meshgrid(mzroi(indMZ),(1:size(Z,2)));
    Z(Z == 0) =nan;

    plot3(X,Y,Z,'LineWidth',1.5)
    ylabel('ATD (bins)')
    xlabel('m/z')
    zlabel('counts')
    if T == 1
    xlim([-.2,.2]+Options.mzTarget)
    end 
%     zL = zlim;
%     zlim([1,max(Z,[],'all')])

    title(['MS ',num2str(T)])
end
% axis tight
%Get location
if Options.getLocation
    datacursormode on
    dcm_obj = datacursormode(fig);
    % Set update function
    set(dcm_obj,'UpdateFcn',@myupdatefcn)
    % Wait while the user to click
    disp('Click line to display a data tip, then press "Return"')

    pause
    % Export cursor to workspace
    info_struct = getCursorInfo(dcm_obj);
    Location(1) = info_struct.Position(2);
    mz =  info_struct.Position(2);

    if isfield(info_struct, 'Position')
        disp('Clicked positioin is')
        disp(info_struct.Position)
    end
else
    Location = nan(1);
end
end
% title(['Mobilogram at ^1t_R: ',num2str(round(info_struct.Position(1),1)),', ^2t_R: ',num2str(round(info_struct.Position(2),1))])%,' at EIC: ',num2str(round(mzroi(a),4))])