combinedWide
res <- pca(combinedWide)


df       <- res$df
scores   <- res$scores
variance <- res$variance



loadings <- as_tibble(reset_stat_defaults()$loadings, rownames = "time") #|> 
  #mutate(time = as.numeric(time))

plot_scores(scores, PC1, PC2, file_path = 'images/pca.png')
plot_loadings(loadings, PC1)
plot_scree(variance)
