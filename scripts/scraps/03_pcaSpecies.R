species_df <- split_by_species(combinedWide)

green <- species_df[['green']]
green <- correct_baseline(green, meta_in = TRUE)
greenRes <- pca(green, center_in = TRUE, scale_in = TRUE)

variance <- greenRes$variance
loadings <- as_tibble(greenRes$loadings, rownames = "time") |> 
  mutate(time = as.numeric(time))
scores <- greenRes$scores

plot_scores(scores, PC1, PC2)
plot_loadings(loadings, PC1)
plot_loadings(loadings, PC2)
plot_scree(variance)




species_df <- split_by_species(combinedWide)

white <- species_df[['white']]

whiteRes <- pca(white)

variance <- whiteRes$variance
loadings <- as_tibble(whiteRes$loadings, rownames = "time") |> 
  mutate(time = as.numeric(time))
scores <- whiteRes$scores

plot_scores(scores, PC1, PC2)
plot_loadings(loadings, PC1)
plot_loadings(loadings, PC2)
plot_scree(variance)
