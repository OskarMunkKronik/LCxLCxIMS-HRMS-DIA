if Options.IMS
    if Options.SaveFig
        % Save figure
        fig = figure(1);
        FigName = [Options.Compound_Name, '_mz',num2str(round(mzroi{k}(a),4)),'Rt_1D',num2str(LocInTime(1)),'_Rt_2D',num2str(LocInTime(2))];
        savefig(fig,[FigName,'_LCxLC.fig'])

        %Save mass spectra vs IMs
        fig = figure(2);
        savefig(fig,[FigName,'_IMSxMZ_MassSpectra.fig'])
        %Save mass spectra
        fig = figure(3);
        savefig(fig,[FigName,'_MassSpectra.fig'])
 
       
        % Save Mass spectra
        mzroi_tmp = mzroi{k};
        save([FigName,'_MassSpectra_tmp_Rt_mzDev.mat'],"MassSpectra_tmp","CompoundInformation","mzroi_tmp")

    end 
end 