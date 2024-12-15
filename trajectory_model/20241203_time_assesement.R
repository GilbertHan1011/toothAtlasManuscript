fullGeneDf <- mes@assays$originalexp@data[1,] %>% as.data.frame()
colnames(fullGeneDf) <- "Gene"
fullGeneDf$Sample <- mes$Project
fullGeneDf$pseudotime <- mes$pseudotime
colnames(fullGeneDf) <- c("y","array_idx","x")
t1 <- Sys.time()
fit_full <- bayesian_gam_regression_nb_shape(x = fullGeneDf$x, y = fullGeneDf$y,array_idx = fullGeneDf$array_idx)
t2 <- Sys.time()
