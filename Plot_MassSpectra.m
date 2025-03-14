function [MassSpectra] = Plot_MassSpectra(mzroi,MSroi,Mod,uRt,idrift,Options,CompoundNumber)
MassSpectra = cell(1,Options.NumTrace);
fig=figure;
    % fig.Position =[2886 209 560 420];
% fig.Position = [-1207 1130  1718 1306];

%
col = 'rb';

for T=1:Options.NumTrace

    if Options.IMS
        StepSize = max(idrift);
        if T==1
            % MassSpecPos= (Options.Location(CompoundNumber,1)-1)*max(uRt)*size(Z,3)+(Options.Location(CompoundNumber,2)-1)*size(Z,3)+Options.Location(CompoundNumber,3) + Options.firstRt;
            MassSpecPos =  find([Mod == Options.Location(CompoundNumber,1) & uRt == Options.Location(CompoundNumber,2) & idrift == Options.Location(CompoundNumber,3)]); %+ Options.firstRt(SampleNumber)
            Sign = 1;
        else
            MassSpecPos= find([Mod == Options.Location(CompoundNumber,1) & uRt == Options.Location(CompoundNumber,2) & idrift == Options.Location(CompoundNumber,3)])+ StepSize;%+ Options.firstRt(SampleNumber);
            Sign = -1;
        end
        
        if Options.IMSCollapse
              if T==1
            % MassSpecPos= (Options.Location(CompoundNumber,1)-1)*max(uRt)*size(Z,3)+(Options.Location(CompoundNumber,2)-1)*size(Z,3)+Options.Location(CompoundNumber,3) + Options.firstRt;
            MassSpecPos =  find([Mod == Options.Location(CompoundNumber,1) & uRt == Options.Location(CompoundNumber,2)]); %+ Options.firstRt(SampleNumber)
            Sign = 1;
        else
            MassSpecPos= find([Mod == Options.Location(CompoundNumber,1) & uRt == Options.Location(CompoundNumber,2)])+ max(idrift);%+ Options.firstRt(SampleNumber);
            Sign = -1;
              end
        end 
    else
        StepSize = 1;
        if T==1
        MassSpecPos     = find([Mod == Options.Location(CompoundNumber,1) & uRt == Options.Location(CompoundNumber,2)],1,"first"); %+ Options.firstRt(SampleNumber);
            Sign = 1;
        else
            MassSpecPos = find([Mod == Options.Location(CompoundNumber,1) & uRt == Options.Location(CompoundNumber,2)],1,"first");% + Options.firstRt(SampleNumber) + 1;
            Sign = -1;
        end
    end
    
    MassSpecPos = MassSpecPos + [-Options.SumNumScan:Options.SumNumScan] * StepSize * Options.NumTrace;
   
    if Options.NormMS
        MassSpectra{T} = sum(MSroi(:,MassSpecPos),2)./max(sum(MSroi(:,MassSpecPos),2));
        ylabel('Counts (normalized)')
    else

        MassSpectra{T} = sum(MSroi(:,MassSpecPos),2);
        ylabel('Counts')
    end

    stem(mzroi,MassSpectra{T}.*Sign,'Marker','none','Color',col(T))
    xlabel('m/z')
    
    hold on
    [~,IndSort]  = sortrows(MassSpectra{T},'descend');

    for n=1:Options.nMZ
        if T == 1
            text(mzroi(IndSort(n)),MassSpectra{T}(IndSort(n))*Sign,num2str(round(mzroi(IndSort(n)),4)),'Rotation',45,'Color',col(T))
        else
            text(mzroi(IndSort(n)),MassSpectra{T}(IndSort(n))*Sign,num2str(round(mzroi(IndSort(n)),3)),'Rotation',315,'Color',col(T))
        end
    end
end
end 