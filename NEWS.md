# fmri.pipeline 0.2-1 (DEVELOPMENT 11Mar2024)

* Fix to use `run_data$tr` inside `build_design_matrix` call
* Always include evt_time 0 in event-aligned outputs for MEDuSA
* Provide `time_audit` argument in event alignment to output additional columns pertaining to alignment timing
* include `volume` as column in outputs of `voxelwise_deconvolution` to support time audits
