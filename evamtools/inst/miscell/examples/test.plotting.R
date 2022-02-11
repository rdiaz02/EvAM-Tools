# library(evamtools)

sample_evam_output <- evam(examples_csd$csd$AND$data)

## Files are named as follows: cpms_orientation_vertex-size_plot-type.png 
png("1_1_all_horizontal_fixed_prob.png", width = 1000, height = 600, units = "px")
plot_CPMs(sample_evam_output, orientation="horizontal", plot_type = "trans_mat", fixed_vertex_size=TRUE)
dev.off()

png("1_2_all_horizontal_obs_prob.png", width = 1000, height = 600, units = "px")
plot_CPMs(sample_evam_output, orientation="horizontal", plot_type = "trans_mat", fixed_vertex_size=FALSE)
dev.off()

png("1_3_all_vertical_obs_prob.png", width = 600, height = 1000, units = "px")
plot_CPMs(sample_evam_output, orientation="vertical", plot_type = "trans_mat", fixed_vertex_size=FALSE)
dev.off()

png("1_4_some_horizontal_obs_prob.png", width = 1000, height = 600, units = "px")
plot_CPMs(sample_evam_output, samples=samples, models=c("OT", "MHN", "CBN"), orientation="horizontal", plot_type = "trans_mat", fixed_vertex_size=FALSE)
dev.off()

png("2_1_all_horizontal_fixed_trm.png", width = 1000, height = 600, units = "px")
plot_CPMs(sample_evam_output, orientation="horizontal", plot_type = "trans_rate_mat", fixed_vertex_size=TRUE)
dev.off()

png("2_2_all_horizontal_obs_trm.png", width = 1000, height = 600, units = "px")
plot_CPMs(sample_evam_output, orientation="horizontal", plot_type = "trans_rate_mat", fixed_vertex_size=FALSE)
dev.off()

png("2_3_all_vertical_obs_trm.png", width = 600, height = 1000, units = "px")
plot_CPMs(sample_evam_output, orientation="vertical", plot_type = "trans_rate_mat", fixed_vertex_size=FALSE)
dev.off()

png("2_4_some_horizontal_obs_trm.png", width = 1000, height = 600, units = "px")
plot_CPMs(sample_evam_output, samples=samples, models=c("OT", "CBN", "CBN"), orientation="horizontal", plot_type = "obs_genotype_transitions", fixed_vertex_size=FALSE)
dev.off()

samples <- evamtools:::sample_CPMs(sample_evam_output, 1000)

png("3_1_all_horizontal_fixed_trans.png", width = 1000, height = 600, units = "px")
plot_CPMs(sample_evam_output, samples=samples, orientation="horizontal", plot_type = "obs_genotype_transitions", fixed_vertex_size=TRUE)
dev.off()

png("3_2_all_horizontal_obs_trans.png", width = 1000, height = 600, units = "px")
plot_CPMs(sample_evam_output, samples=samples, orientation="horizontal", plot_type = "obs_genotype_transitions", fixed_vertex_size=FALSE)
dev.off()

png("3_3_all_vertical_obs_trans.png", width = 600, height = 1000, units = "px")
plot_CPMs(sample_evam_output, samples=samples, orientation="vertical", plot_type = "obs_genotype_transitions", fixed_vertex_size=FALSE)
dev.off()

png("3_4_some_horizontal_obs_trans.png", width = 1000, height = 600, units = "px")
plot_CPMs(sample_evam_output, samples=samples, models=c("OT", "CBN"), orientation="horizontal", plot_type = "obs_genotype_transitions", fixed_vertex_size=FALSE)
dev.off()