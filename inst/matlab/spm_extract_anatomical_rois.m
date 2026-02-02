function spm_extract_anatomical_rois(cfg)
    % Extracts ROIs from one subject.
    % Arguments:
    %       cfg: A struct with the following fields
    %
    %       target_dir: path to subject's first level analysis
    %                   results.
    %
    %       adjust:     0 to not adjust, NaN to adjust
    %                   for everything, else index of F-contrast to use.
    %
    %       contrast:   Index of contrast of interest.
    %
    %       threshold:  The p value at which to threshold the
    %                   contrast image.
    %
    %       threshdesc: Whether to correct for  comparisons.
    %                   'none' for no correction, 'FWE' for correction.
    %
    %       extent:     Exclude clusters with number of significant
    %                   voxels that is less than extent.
    %
    %       name:       Name of VOI output.
    %                   'VOI_' is prefixed to it by spm_regions.m
    %
    %       mask:       String paths to ROI mask.
    %

    spm('defaults', 'FMRI');
    spm_jobman('initcfg');

    % spm_regions needs to be in SPM.mat directory.
    starting_dir = pwd;
    cd(cfg(1).target_dir);

    for ii = 1 : numel(cfg)
        jobs{ii}.spm.util.voi.spmmat = {fullfile(...
            cfg(1).target_dir, 'SPM.mat')};
        jobs{ii}.spm.util.voi.adjust = cfg(ii).adjust;
        jobs{ii}.spm.util.voi.session = cfg(ii).session;
        jobs{ii}.spm.util.voi.name = cfg(ii).name;

	%use contrast for defining region
	if isfield(cfg(ii), 'contrast') && ~isempty(cfg(ii).contrast)
	  % The image data
          jobs{ii}.spm.util.voi.roi{1}.spm.spmmat = {''};
          jobs{ii}.spm.util.voi.roi{1}.spm.contrast = cfg(ii).contrast;
          jobs{ii}.spm.util.voi.roi{1}.spm.threshdesc = cfg(ii).threshdesc;
          jobs{ii}.spm.util.voi.roi{1}.spm.thresh = cfg(ii).threshold;
          jobs{ii}.spm.util.voi.roi{1}.spm.extent = cfg(ii).extent;

	  % The anatomical ROI mask
          jobs{ii}.spm.util.voi.roi{2}.mask.image = {cfg(ii).mask};

	  % The logical expression combining images.
          jobs{ii}.spm.util.voi.expression = 'i1 & i2';
	else
	  % Use only an anatomical ROI mask
          jobs{ii}.spm.util.voi.roi{1}.mask.image = {cfg(ii).mask};

	  % The logical expression combining images.
          jobs{ii}.spm.util.voi.expression = 'i1';
	end
	
    end

    try
        spm_jobman('run', jobs);
    catch ME
        cd(starting_dir);
        rethrow(ME);
    end
    cd(starting_dir)
end
