%% S_PlotLCxLCxIMS-HRMS
%% Load data files (ONLY RUN THIS WHEN YOU OPEN MATLAB)

% Initialize
OptionsImport.Root_dir = 'C:\Users\mht541\Documents\2D_LC_IMS_R\2DIMS'; % Change the first part to make the folder '2DIMS' the root directory
OptionsImport.Save_fig_folder = 'C:\Users\mht541\Documents\2D_LC_IMS_R\2DIMS\Figures'
cd(OptionsImport.Root_dir)
ProjectPBA
%%

% Initialize
OptionsImport.Root_dir = 'C:\Users\mht541\Documents\2D_LC_IMS_R\2DIMS\'; % Change the first part to make the folder '2DIMS' the root directory
cd(OptionsImport.Root_dir)
ProjectPBA
OptionsImport.Save_fig_folder='C:\Users\mht541\Documents\2D_LC_IMS_R\SuspectScreening'; %Change this to the folder where the figures should be saved

%Sample names
%RPLCxHILIC
SampleName = {'RPLCxHILIC - Neg','RPLCxHILIC - Pos','RPLCxHILICxcIMS - Neg','RPLCxHILICxcIMS - Neg - Standards'};% Write in your sample names

% Options for plotting
Options = struct('IMS',true, ...    true if IMS 
    'IMSCollapse',false, ...        
    'NormMS',false, ...              %true if mass spectra should be normalized to the base peak in the MS1 (low energy) and MS2 (high energy) trace respectively
    'PlotType','EIC', ...           %Options: 1)Total ion chromatogram: TIC, 2) base peak chromatogram: BPC, 3) extracted ion chromatogram (EIC)
    'NumTrace',2, ...               % If you have MSe data this should be 2, if only MS write 1
    'nMZ',10, ...                   %The number of m/z that you want to annotate in the mass spectra in the the MS1 (low energy) and MS2 (high energy) trace respectively
    'Weight',1, ...                 %This scales the intensities in the chromatogram to '^1/Weight' to better see low intensities 
    'mzTarget', 457.077635, ...       %Change this to the correct m/z 
    'Compound_Name','quercetin',... %Change this to save the figure under the correct folder
    'PlotIMS',false, ...            
    'PlotMS',false, ...
    'Plot2D','surf', ...
    'fillIn',nan, ...,              
    'SumNumScan',3,...              %+/-The number of scan points which should be used to sum into the mass spectra
    'getLocation',true, ...         %Should be true if you choose the position to take the LCxIMS and mass spectra from.
    'getLocationAuto', false, ...
    'SaveFig',true);                % true: If you want to save figures of the chromatogram, LCxIMS, and mass spectra (.fig) and the mass spectra as a .mat file 

%%
%Don't change! 
modTimeVec = repelem(24,length(fileList)); %Write in the modulation time
PhaseShift = [0,0,0,0]; %A positive number is how much you will shift the phase

%Sample to plot corresponds to the number in SampleName: 
SampleVec = 3; %sample to plot
Options.Compound_Name = 'test'; %write in your compound name 
Options.mzTarget =301.0345;% Write in the exact m/z 
% Options_ID.targetInChIKey = CompoundsInformation{n,2}{:};
% 


for k =[SampleVec]

    modTime= modTimeVec(k);
    Options.IMS = IMSVec(k);
    if ~Options.IMS 
        Dt{k} = ones(size(MSroi{k},2),1);
%         Options_Align.mzTarget = 157.9284;%297.839;
    end 
    if exist([OptionsImport.Save_fig_folder,'\DataFiles\Figures\',fileList(k).name(1:end-4),'\',Options.Compound_Name])
        cd([OptionsImport.Save_fig_folder,'\DataFiles\Figures\',fileList(k).name(1:end-4),'\',Options.Compound_Name])
    else
        mkdir([OptionsImport.Save_fig_folder,'\DataFiles\Figures\',fileList(k).name(1:end-4),'\',Options.Compound_Name])
        cd([OptionsImport.Save_fig_folder,'\DataFiles\Figures\',fileList(k).name(1:end-4),'\',Options.Compound_Name])
  
    end

    % Adjust the phase shift.
%     if k == 1
%         Options.firstRt(k) =  find(Rt{k} >= PhaseShift(k)  & trace{k} == 1,1,'first');
%     end 
    Options.firstRt(k)      =   find(Rt{k} >= PhaseShift(k)  & trace{k} == 1,1,'first');
    
    % Adjust the Rt vector
    Rt2{k} = Rt{k}(Options.firstRt(k):end);
    Rt2{k} = Rt2{k}-Rt2{k}(1);
    

    %Find modulateoions and 2D Rts
    [Mod,uRt] = FoldChrom(round(Rt2{k},2),modTime,trace{k}(Options.firstRt(k):end));

    %Drift times
    [udrift,ia,idrift] = unique(Dt{k}(Options.firstRt(k):length(Dt{k})));

      % Plotting
    S = 1;
    close all
    [minMZ,a] = min(abs(mzroi{k}-Options.mzTarget));
%     a = a-1;
    if  (minMZ/Options.mzTarget*10^6) <= 1000

        %Print sample name and mass accuracy
        fprintf(['Sample name: ',SampleName{k},'\n'])
        fprintf(['Compound name: ',Options.Compound_Name,'\n'])
        ppmDev = (mzroi{k}(a)-Options.mzTarget)/Options.mzTarget*10^6;
        fprintf(['Target m/z: ',num2str(Options.mzTarget),', ppm deviation: ',num2str(ppmDev),'\n'])
        
        if Options.IMS
        
            Ind = sum(MSroi{k}(:,1:find(Rt{k}>modTimeVec(k)*10,1,'first'))>50,2)>0; %
        else 

            Ind = sum(MSroi{k}(:,1:find(Rt{k}>modTimeVec(k)*10,1,'first'))>500,2)>0; %
        end 

        [X,Y,Z,Options.Location(S,1:2),LocInTime(1:2),indTrace,Intensity] = PlotLCxLC(mzroi{k},MSroi{k}(:,Options.firstRt(k):end),Mod,Rt2{k},uRt,idrift,modTime,Options,trace{k}(Options.firstRt(k):end),[]);
        
        title(['EIC, mz:  ',num2str(round(Options.mzTarget,4))])
        nn = 1
%         [width,roots] = PeakWidth(timeVector,intensityVector,Apex,percentHeight)
        if Options.IMS
            % IMS
            [Options.Location(S,3)] = Plot_IMSxMS(mzroi{k},MSroi{k}(:,Options.firstRt(k):end),Mod,uRt,Options,S,trace{k}(Options.firstRt(k):end));
        end

        % Mass Accuracy
        mzAccuracy(S) = ppmDev;
        [MassSpectra_tmp] = Plot_MassSpectra(mzroi{k},MSroi{k}(:,Options.firstRt(k):end),Mod,uRt,idrift,Options,S);
        [MassSpectra_tmp_clean] = Plot_MassSpectra(mzroi{k}(~Ind),MSroi{k}(~Ind,Options.firstRt(k):end),Mod,uRt,idrift,Options,S);
   
        MassSpectra(S,:) = MassSpectra_tmp;
        MassSpectra_Cleaned(S,:) = MassSpectra_tmp_clean;

        if Options.IMS
            CompoundInformation =  [LocInTime,Options.Location(S,3),mzroi{k}(a),mzAccuracy(S)]
        else
            [LocInTime,mzroi{k}(a),mzAccuracy(S)]
        end

        pause
        %Save Figures   
        SaveFigures_2DLC_IMS

    end
end
