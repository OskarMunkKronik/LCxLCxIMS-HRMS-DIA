function [X,Y,Z,Location,LocInTime,indTrace,intensity] = PlotLCxLC(mzroi,MSroi,Mod,timeVector,uRt,idrift,modTime,Options,trace,Filter)
% close all
if isequal(Options.PlotType,'TIC')
    intensityVector = sum(MSroi,1);
end
if isequal(Options.PlotType,'BPC')
    intensityVector = max(MSroi,[],1);
end
if isequal(Options.PlotType,'EIC')
    [~,a] = min(abs(mzroi-Options.mzTarget));
    intensityVector = sum(MSroi(a-1:a+1,:),1);
end
Z=accumarray([Mod,uRt,idrift,trace],full(intensityVector),[max(Mod),max(uRt),max(idrift),max(trace)],[],Options.fillIn);
fig = figure;
% fig.Position = [-3119 1240 1902 1057];
indTrace=unique(uRt(trace == 1));

Rt1D = [0:size(Z,1)-1].*modTime/60;
Rt2D=[1:max(indTrace)].*modTime./max(indTrace);
% [Y,X]   = meshgrid(Rt2D(indTrace),Rt1D);
[Y,X]   = meshgrid(Rt2D,Rt1D);
surf(X(:,indTrace),Y(:,indTrace),squeeze(sum(Z(:,indTrace,:,1),3,'omitnan')).^(1/Options.Weight),'EdgeColor','interp','FaceColor','interp')
% xlim([8 60])
% Aestethics
xlabel('^1t_R (min)')
ylabel('^2t_R (sec)')
colorbar
view([0 90])
axis tight
intensity =[];


%Get location
if Options.getLocation
    if Options.getLocationAuto
        [Location]  =   FastPeakFind(squeeze(sum(Z(:,indTrace,:,1),3))', 100, Filter);
        Location    =   reshape(Location,2,[])';

        info_struct.Position = zeros(size(Location));
        info_struct.Position(:,1) = (Location(:,1)-1).*modTime/60;
        Location(:,2) = indTrace(Location(:,2));
        for n = 1:size(Location,1)
        info_struct.Position(n,2) = Y(Location(n,1),Location(n,2));
        end 
    else

        datacursormode on
        dcm_obj = datacursormode(fig);
        % Set update function
        set(dcm_obj,'UpdateFcn',@myupdatefcn)
        % Wait while the user to click
        disp('Click line to display a data tip, then press "Return"')

        pause
        % Export cursor to workspace
        info_struct = getCursorInfo(dcm_obj);
    end
    for n = 1:size(info_struct.Position,1)
        Location(n,1)=find(Rt1D == info_struct.Position(n,1));
        Ind2D =  find(Y(Location(n,1),:) == info_struct.Position(n,2));
        
    
    % if length(Ind2D) > 1
    [~,a] = max(sum(Z(Location(n,1),Ind2D,:,1),3));
    Location(n,2) = Ind2D(a);
    % = ; %find(Rt2D == info_struct.Position(2));
    LocInTime(n,:) = info_struct.Position(n,1:2);
    intensity = info_struct.Position(n,3);
    end
    
    if isfield(info_struct, 'Position')
        disp('Clicked positioin is')
        disp(info_struct.Position)
    end
else
    [LocInTime,Location] = deal(nan(1,2));
end


end