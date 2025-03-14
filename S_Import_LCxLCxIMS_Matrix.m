%% Addpaths
OptionsImport.BaseDir = 'C:\Users\mht541\Documents\2D_LC_IMS_R\2DIMS'; %Write the path to your folder where you have the data
OptionsImport.CDF_dir = 'F:\PhD\Matrix\IMS\RPLCxHILIC_CDF\Matrix\RPLCxHILIC\SecondAttempt'; %write the path to the cdf files  
OptionsImport.Mat_dir = 'C:\Users\mht541\Documents\2D_LC_IMS_R\2DIMS\Matrix\DataFiles\Mat\SecondAttempt'; %Write the path to where you want the processed data to be stored
SampleSetSaveName = '250227_RPLCxHILIC_Negative_NoFilter_NoIMS_test'; % Write in the sampleset save name
cd(OptionsImport.BaseDir)
ProjectPBA

% Options files
OptionsROI = struct('mzerror',14,...
    'ppm',true,...
    'thresh',10,...
    'wmean',true,...
    'GapAllowed',2,...
    'minroi',4,...
    'NumTrace',2,...
    'verbose',true,...
    'prefilter',true,...
    'fillIn',0, ...
    'IMS',true, ...
    'CollapseRoIs', true,...
    'RefMass',true,...
    'RefMass_mz',    554.2620,...
    'OneFile',true,...
    'RefMass_dir','F:\PhD\Matrix\IMS\RPLCxHILIC_CDF\Matrix\RPLCxHILIC\SecondAttempt\MS3');

%LockMass
cd(OptionsImport.CDF_dir)
fileList = dir('*CDF');
[mzroi, MSroi, RawData,runTime,sizeMZRoI,Rt,Dt,trace] = deal(cell(length(fileList),1)); %Allocates cells
SampleVec = [1:length(fileList)]; % Selects the samples to process 
SampleVec_mz_determination = [3] % Selects the sample to use for m/z deviation optimization
IMSVec = repelem(true,length(fileList)); %true means that it is an IMS file
SampleName = {'RPLCxHILIC - Neg','RPLCxHILIC - Pos','RPLCxHILICxcIMS - Neg','RPLCxHILICxcIMS - Neg - Standards'};%'RPLCxHILICxcIMS - Pos',


%% Optimize mzerror
clf
  MassAccuracy = 2; % The lowest m/z deviation test 
    MZmultiplyFactors = [1:2:8,9,10]; % The multiples tested for the the mass accuracy 
 
    %Optimization of the m/z deviation
mzDev = nan(length(SampleVec_mz_determination),1)
for n = 1:length(SampleVec_mz_determination)
    k = SampleVec_mz_determination(n);
    fprintf([fileList(k).name,'\n'])

     OptionsROI.IMS
    [mzDev(n),minroi,nmz]=OptParamRoi(fileList(k).name,MassAccuracy,MZmultiplyFactors,OptionsROI.minroi,OptionsROI);
    % pause
end
OptionsROI.mzerror = round(mean(mzDev,'omitnan'));
%%
% Make fileList

% Performs ROI detection of the samples
for k =SampleVec
    cd(OptionsImport.CDF_dir)
    fprintf([fileList(k).name,'\n'])
    OptionsROI.IMS = IMSVec(k);
 
    if IMSVec(k)
        [mzroi{k},MSroi{k},~,         runTime{k},sizeMZRoI{k},Rt{k},Dt{k},trace{k}]=ROIpeaks_ACN(fileList(k).name,OptionsROI);%thresh,mzerror,minroi,sScan,eScan,ppm,wmean);
        
        nBins = length(unique(Dt{k}));
        max_scan = size(MSroi{k},2);
        
        %Lower to max_scan
        MSroi{k} = MSroi{k}(:,1:max_scan);
        Rt{k}    = Rt{k}(1:max_scan);
        Dt{k}    = Dt{k}(1:max_scan);
        trace{k} = trace{k}(1:max_scan)';
        firstRt(k) = find(cumsum(sum(MSroi{k}))'>0 & trace{k} == 1 ,1,'first');
        Dt{k} = Dt{k}(firstRt(k):end);
    else
    
        [mzroi{k},MSroi{k},~,runTime{k},sizeMZRoI{k},Rt{k},~,trace{k}]=ROIpeaks_ACN(fileList(k).name,OptionsROI);%thresh,mzerror,minroi,sScan,eScan,ppm,wmean);


        %Lower to max_scan
        max_scan = floor(size(MSroi{k},2)/2)*2;
        MSroi{k} = MSroi{k}(:,1:max_scan);
        Rt{k}    = Rt{k}(1:max_scan);
        trace{k} = trace{k}(1:max_scan)';
        firstRt(k) = find(sum(MSroi{k})'> 1000 & trace{k} == 1 ,1,'first');%find(cumsum(sum(MSroi{k}))'>0 & Dt{k} == 0 & trace{k} == 1 ,1,'first');
      
    end
    
    trace{k} =trace{k}(firstRt(k):end);
    Rt{k} = Rt{k}(firstRt(k):end);
    Rt{k} = Rt{k}-Rt{k}(1);

    MSroi{k}=MSroi{k}(:,firstRt(k):end);
end

%Saves the data 
cd(OptionsImport.Mat_dir)
save([SampleSetSaveName,'.mat'],'-v7.3')
