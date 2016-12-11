library(normalmultinomial)

if(!exists('build')) build = 'N_01000-n_00050-s_00001-seed_00001'
pattern_build = "N_([0-9]+)-n_([0-9]+)-s_([0-9]+)-seed_([0-9]+)*"

load(sprintf('datasets/dataset-%s.RData', build))

fit <- dm_fit(XZ)

save(fit, file = sprintf('datasets/replacement-%s-method_dm.RData', build))
