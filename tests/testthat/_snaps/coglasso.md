# print.coglasso works

    Code
      cg <- coglasso(multi_omics_sd_micro, p = 4, nlambda_w = 3, nlambda_b = 3, nc = 2,
        verbose = FALSE)
      print(cg)
    Output
      Networks estimated with collaborative graphical lasso
      
      The call was:
      coglasso(data = multi_omics_sd_micro, p = 4, nlambda_w = 3, nlambda_b = 3, 
          nc = 2, verbose = FALSE)
      
      The total number of hyperparameter combinations explored was:
      18
      The values explored for lambda within were:
      0.8743, 0.2765, 0.0874
      The values explored for lambda between were:
      0.4816, 0.1523, 0.0482
      The values explored for c were:
      1, 10
      
      Networks are made of 2 omics layers, for a total of 6 nodes
      For each layer they have: 4 and 2 nodes, respectively
      
      Select the best network with:
      sel_cg <- select_coglasso(cg)

