# print.select_coglasso works with all functions

    Code
      cg <- coglasso(multi_omics_sd_micro, p = 4, nlambda_w = 3, nlambda_b = 3, nc = 2,
        verbose = FALSE)
      sel_cg <- select_coglasso(cg, method = "ebic", verbose = FALSE)
      print(sel_cg)
    Output
      Selected network estimated with collaborative graphical lasso
      
      The call was:
      select_coglasso(coglasso_obj = cg, method = "ebic", verbose = FALSE)
      
      The model selection method was:
      ebic
      The density of the selected network is:
      0.4666667
      
      Networks are made of 2 omics layers, for a total of 6 nodes
      For each layer they have: 4 and 2 nodes, respectively
      
      The selected value for lambda within is:
      0.0874
      The selected value for lambda between is:
      0.4816
      The selected value for c is:
      0.01
      
      The total number of hyperparameter combinations explored was:
      18
      The values explored for lambda within were:
      0.8743, 0.6852, 0.0874
      The values explored for lambda between were:
      0.4816, 0.3775, 0.0482
      The values explored for c were:
      0.01, 100
      
      Plot the selected network with:
      plot(sel_cg)
    Code
      sel_cg <- select_coglasso(cg, method = "xestars", rep_num = 3, verbose = FALSE)
      print(sel_cg)
    Output
      Selected network estimated with collaborative graphical lasso
      
      The call was:
      select_coglasso(coglasso_obj = cg, method = "xestars", rep_num = 3, 
          verbose = FALSE)
      
      The model selection method was:
      xestars
      The density of the selected network is:
      0
      
      Networks are made of 2 omics layers, for a total of 6 nodes
      For each layer they have: 4 and 2 nodes, respectively
      
      The selected value for lambda within is:
      0.8743
      The selected value for lambda between is:
      0.4816
      The selected value for c is:
      0.01
      
      The total number of hyperparameter combinations explored was:
      18
      The values explored for lambda within were:
      0.8743, 0.6852, 0.0874
      The values explored for lambda between were:
      0.4816, 0.3775, 0.0482
      The values explored for c were:
      0.01, 100
      
      Plot the selected network with:
      plot(sel_cg)
    Code
      sel_cg <- xstars(cg, rep_num = 3, verbose = FALSE)
      print(sel_cg)
    Output
      Selected network estimated with collaborative graphical lasso
      
      The call was:
      xstars(coglasso_obj = cg, rep_num = 3, verbose = FALSE)
      
      The model selection method was:
      xstars
      The density of the selected network is:
      0
      
      Networks are made of 2 omics layers, for a total of 6 nodes
      For each layer they have: 4 and 2 nodes, respectively
      
      The selected value for lambda within is:
      0.8743
      The selected value for lambda between is:
      0.4816
      The selected value for c is:
      0.01
      
      The total number of hyperparameter combinations explored was:
      18
      The values explored for lambda within were:
      0.8743, 0.6852, 0.0874
      The values explored for lambda between were:
      0.4816, 0.3775, 0.0482
      The values explored for c were:
      0.01, 100
      
      Plot the selected network with:
      plot(sel_cg)
    Code
      sel_cg <- bs(multi_omics_sd_micro, p = 4, nlambda_w = 3, nlambda_b = 3, nc = 2,
        rep_num = 3, verbose = FALSE)
      print(sel_cg)
    Output
      Selected network estimated with collaborative graphical lasso
      
      The call was:
      bs(data = multi_omics_sd_micro, p = 4, nlambda_w = 3, nlambda_b = 3, 
          nc = 2, rep_num = 3, verbose = FALSE)
      
      The model selection method was:
      xestars
      The density of the selected network is:
      0
      
      Networks are made of 2 omics layers, for a total of 6 nodes
      For each layer they have: 4 and 2 nodes, respectively
      
      The selected value for lambda within is:
      0.8743
      The selected value for lambda between is:
      0.4816
      The selected value for c is:
      0.01
      
      The total number of hyperparameter combinations explored was:
      18
      The values explored for lambda within were:
      0.8743, 0.6852, 0.0874
      The values explored for lambda between were:
      0.4816, 0.3775, 0.0482
      The values explored for c were:
      0.01, 100
      
      Plot the selected network with:
      plot(sel_cg)

