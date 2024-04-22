# rm(list = ls())
library(magrittr)
#### Functions ####


add_dependencies = function(Balancing_Method) {
  # for each nuisance function of each method, add the required number of split
  # and nuisance functions needed to be called to compute the nuisance function
  
  lapply(
    Balancing_Method,
    function(met) {
      # for each nuisance function get the names of the other nuisance function
      # required to compute the nuisance function
      args_name = lapply(
        met[["nuisance_fct"]],
        function(nui) {
          formalArgs( nui[["train"]] ) %>%
            `[`( . %in% names(met[["nuisance_fct"]]) ) %>%
            { if (length(.) == 0) NULL else . }
        })
      
      # if a nuisance function do not require other nuisance function then it
      # only requires one split and it does not have any dependencies
      nbr_required_split = sapply(args_name, function(name) if (is.null(name)) 1 else NA)
      dependencies = lapply(args_name, function(name) if (is.null(name)) list() else NA)
      
      
      while (any(is.na(nbr_required_split))) {
        # get the names of nuisance function  for which 'nbr_required_split' and
        # 'dependencies' are still not computed
        to_compute = nbr_required_split %>%
          { names(.)[is.na(.)] }
        
        for (name in to_compute) {
          # the number of split required by a nuisance function is the sum of
          # required number split of its direct dependencies + 1
          nbr_required_split[name] = args_name[[name]] %>%
            nbr_required_split[.] %>%
            sum(1)
          
          # the dependencies of a nuisance function is the union of the
          # dependencies's of the nuisance function it requires to be computed
          dependencies[[name]] = args_name[[name]] %>%
            met[["nuisance_fct"]][.] %>%
            lapply(`[[`, "train") %>%
            c(dependencies[[args_name[[name]]]]) %>%
            unique() %>%
            { if (any(is.na(.))) NA else . }
        }
      }
      
      # save the values
      for (nui_name in names(met[["nuisance_fct"]])) {
        met[["nuisance_fct"]][[nui_name]]["nbr required split"] = nbr_required_split[nui_name]
        met[["nuisance_fct"]][[nui_name]][["dependencies"]] = dependencies[[nui_name]]
      }
      return(met)
    }
  )
}

add_unique_nuisance = function(Balancing_Method) {
  # create the skeleton for 'unqiue_nuisance', a list containing the minimum
  # number of nuisance function needed to be computed to crossfit the methods in
  # 'Balancing_Method'.
  #
  # for each nuisance function, 'unique_nuisance' contains :
  #   - the nuisance function;
  #   - the name given to this nuisance function by the 1st method in Balancing_Method
  #     using the nuisance function;
  #   - the names of the other nuisance functions required to compute the nuisance
  #     function (the names are also retrived from the 1st method calling the nuisance function),
  #   - the index in 'unique_nuisance' of the other nuisance functions required to compute
  #     the nuisance function;
  #   - a boolean equal to TRUE if the nuisance function should predict the values in
  #     the other splits;
  #   - a boolean equal to TRUE if the nuisance function should predict the values in test_dataset;
  #   - the indexes of the methods requiring the nuisance function.
  
  # check if the number of split for each method is the same
  sapply(Balancing_Method, `[[`, "nbr_split") %>%
    { any(.[1] != .) } %>%
    if (.) stop("'add_unique_nuisance' should only be used for 'Balancing_Method'
                containing method with the same 'nbr_split'")
  
  nbr_split = Balancing_Method[[1]][["nbr_split"]]
  
  # add dependencies if missing
  Balancing_Method = sapply(
    Balancing_Method,
    function(met) {
      sapply(met[["nuisance_fct"]], `[[`, "nbr required split") %>%
        sapply(is.null) %>% any()
    }) %>% any() %>%
    if (.) add_dependencies(Balancing_Method) else Balancing_Method
  
  # list of unique nuisance function order by the number of split required
  unique_nuisance = lapply(Balancing_Method, `[[`, "nuisance_fct") %>%
    # concatenate the nuisance functions
    do.call(c,.) %>% 
    # remove duplicate
    unique() %>% 
    # order the list
    `[`(sort.list(sapply(., `[[`, "nbr required split")))
  
  
  Balancing_Method[["unique_nuisance"]] = lapply(
    unique_nuisance,
    function(nui) {
      # get the name of 'nui' and a method calling 'nui'. 
      # get the index of the method requiring 'nui'
      met_idx = sapply(
        seq_along(Balancing_Method),
        function(i) {
          match(list(nui), Balancing_Method[[i]][["nuisance_fct"]]) %>%
            sapply(is.na) %>% all(.)
        }) %>% `!` %>% which()
      met = Balancing_Method[[met_idx[1]]]
      nui_name = match(list(nui), met[["nuisance_fct"]]) %>%
        `[`(names(met[["nuisance_fct"]]), .) %>% .[1]
      
      
      # should 'nui' predict the value of the test dataset
      predict_test_dataset = sapply(
        Balancing_Method,
        function(met) {
          if (!(nui_name %in% formalArgs(met[["method"]]))) return(FALSE)
          return(list(nui) %in% met[["nuisance_fct"]])
        }) %>% any()
      
      
      # should 'nui' predict on the other splits
      predict_split = sapply(
        unique_nuisance,
        function(other_nui) {
          if (!(nui_name %in% formalArgs(other_nui[["train"]]))) return(FALSE)
          return(nui["train"] %in% other_nui[["dependencies"]])
        }) %>% any()
      
      
      # names of the required nuisance functions for 'nui' 
      names_required_nui = names( met[["nuisance_fct"]] ) %>%
        `[`( . %in% formalArgs(nui[["train"]]) ) %>%
        { if (length(.) == 0) NULL else . }
      
      # idx of the required nuisance function in unique_nuisance
      nui_unique_nuisance_idx = if (is.null(names_required_nui)) NULL else sapply(
        seq_along(names_required_nui),
        function(j) {
          # get the name of the required nuisance function
          met[["nuisance_fct"]][names_required_nui][[j]][c("train", "dependencies")] %>%
            list() %>%
            match(lapply(unique_nuisance, `[`, c("train", "dependencies")))
        })
      
      
      lapply(
        1:nbr_split,
        function(i) {
          
          if ( length(nui[["dependencies"]]) == 0 ) return(NULL)
          
          # get the index of the split used to train the required nuisance
          # functions
          names_required_nui %>%
            { met[["nuisance_fct"]][.] } %>%
            sapply(`[[`, "nbr required split") %>%
            c(0,.) %>% cumsum() %>% .[-length(.)] %>%
            { i - . - 1 } %>%
            { (. - 1) %% nbr_split + 1} %>%
            list("nui_train_split_idx" = .)
        }) %>%
        c(nui,
          list("nui_name" = nui_name,
               "ignore" = FALSE,
               "met_idx" = met_idx,
               "predict_split" = predict_split,
               "predict_test_dataset" = predict_test_dataset,
               "names_required_nui" = names_required_nui,
               "nui_unique_nuisance_idx" = nui_unique_nuisance_idx)
        )
    })
  return(Balancing_Method)
}


add_method_info = function(Balancing_Method) {
  # creat 'method_info', a list containing useful information (such as the index
  # in 'unique_nuisance' of the nuisance function needed to compute a method) to
  # crossfit a method in 'Balancing_Method'.
  
  # check if 'unique_nuisance' is already in 'Balancing_Method'
  if (is.null(Balancing_Method[["unique_nuisance"]]))
    Balancing_Method %<>% add_unique_nuisance()
  
  unique_nuisance = Balancing_Method[["unique_nuisance"]]
  
  nbr_split = Balancing_Method[[1]][["nbr_split"]]
  Balancing_Method[["method_info"]] = lapply(
    Balancing_Method[-length(Balancing_Method)],
    function(met) {
      
      # names of the nuisance functions needed to be passed on as an arg to
      # met[["method"]]
      names_required_nui = names(met[["nuisance_fct"]]) %>%
        `[`( . %in% formalArgs(met[["method"]]) ) %>%
        { if (length(.) == 0) NULL else . }
      
      # index in 'unique_nuisance' of the nuisance functions needed to be passed
      # on as an arg to met[["method"]]
      nui_unique_nuisance_idx = if (is.null(names_required_nui)) NULL else sapply(
        seq_along(names_required_nui), 
        function(i) {
          met[["nuisance_fct"]][names_required_nui][[i]][c("train", "dependencies")] %>%
            list() %>%
            match( lapply(unique_nuisance,`[`, c("train", "dependencies")) )
        })
      
      
      lapply(
        1:nbr_split,
        function(split_idx) {
          # if met does
          if (is.null(names_required_nui)) return(NULL)
          
          # get the index of the split used to train the required nuisance
          # functions
          met[["nuisance_fct"]][names_required_nui] %>%
            sapply(`[[`, "nbr required split") %>%
            c(0, .) %>% cumsum() %>% `[`( -length(.) ) %>%
            { split_idx - . - 1 } %>%
            { (. - 1) %% nbr_split + 1} %>%
            list("nui_train_split_idx" = .)
        }) %>%
        c(list("names_required_nui" = names_required_nui,
               "nui_unique_nuisance_idx" = nui_unique_nuisance_idx)
        )
    })
  return(Balancing_Method)
}


subset_Balancing_Method = function(Balancing_Method, idx) {
  to_ignore = sapply(
    lapply(Balancing_Method[["unique_nuisance"]], `[[`, "met_idx"),
    function(i) !(any(i %in% idx))
  ) %>% which()
  
  for (i in to_ignore) Balancing_Method[["unique_nuisance"]][[i]][["ignore"]] = TRUE
  
  c(Balancing_Method[idx],
    list("method_info" = Balancing_Method[["method_info"]][idx],
         "unique_nuisance" = Balancing_Method[["unique_nuisance"]])
  )
}


simple_crossfit = function(X_tr, Y_tr, Z_tr, X_te, Y_te, Z_te,
                           Balancing_Method, nbr_split) {
  n = nrow(X_tr) # number of observations in train dataset
  p = ncol(X_tr) # number of covariates
  fail_output = structure("", class = "try-error")
  
  if (is.null(Balancing_Method[["method_info"]]))
    Balancing_Method %<>% add_method_info
  
  unique_nuisance = Balancing_Method[["unique_nuisance"]]
  method_info = Balancing_Method[["method_info"]]
  Balancing_Method = Balancing_Method[-length(Balancing_Method) + 0:1]
  
  test_dataset = list(
    "X" = X_te,
    "Y" = Y_te,
    "Z" = Z_te)
  
  # split the data evenly without overlap
  splited_data = data.frame(X_tr, Y_tr, Z_tr) %>%
    # shuffle the rows
    .[sample(n), ] %>%
    # split [X_tr, Y_tr, Z_tr] into "nbr_split" splits
    split( rep_len(1:nbr_split, n) ) %>%
    # unbind [X_tr, Y_tr, Z_tr]
    lapply( function(data) list(
      "X" = data[, 1:p],
      "Y" = data[, p + 1],
      "Z" = data[, p + 2])
    )
  
  # train the models on each split and predict the values on the other splits
  for (nui_idx in seq_along(unique_nuisance)) {
    
    nui = unique_nuisance[[nui_idx]]
    
    # skip nuisance function not needed by the methods
    if (nui[["ignore"]]) next
    
    unique_nuisance[[nui_idx]] = lapply(
      1:nbr_split,
      function(i) {
        
        out = nui[i]
        
        # get the predicted values on the split i of the required nuisance
        # functions trained on the split nui_train_split_idx
        required_nui_input =  lapply(
          seq_along(nui[["names_required_nui"]]),
          function(j) {
            # get the prediction on the split i
            unique_nuisance %>%
              `[[`(nui[["nui_unique_nuisance_idx"]][j]) %>%
              `[[`(nui[[i]][["nui_train_split_idx"]][j]) %>%
              `[[`("predicted.value") %>%
              `[[`(i)
          }) %>%
          `names<-`( nui[["names_required_nui"]] ) %>%
          do.call(c, .)
        
        
        # train the nuisance function on the split i
        out[["trained.model"]] = splited_data[[i]] %>%
          c(required_nui_input) %>%
          # remove unused args
          `[`( names(.) %in% formalArgs(nui[["train"]]) ) %>%
          { try( do.call(nui[["train"]], .) ) }
        if (length(out[["trained.model"]])==0) browser()
        
        out[["fail"]] = inherits(out[["trained.model"]], "try-error")
        
        # predicted the values on the other splits and test_dataset
        out[["predicted.value"]] = lapply(
          1:(nbr_split + 1),
          function(j) {
            # do not try to predict if the training failed
            if (out[["fail"]]) return(out[["trained.model"]])
            # if there is no need to predict on the splits
            if (j <= nbr_split & !nui[["predict_split"]]) return(NULL)
            # if there is no need to predict on the test dataset
            if (j > nbr_split & !nui[["predict_test_dataset"]]) return(NULL)
            # do not predict on the split used to train the model
            if (j==i) return(NULL)
            
            # predict on the other splits and the testing data set
            (if (j == nbr_split+1) test_dataset else splited_data[[j]]) %>%
              c(out["trained.model"]) %>%
              # remove unused args
              .[names(.) %in% formalArgs(nui[["predict"]])] %>%
              # preidict the values in the current split
              { try(do.call(nui[["predict"]], .)) }
          }) %>%
          # rename the lists containing the predictions to 1, 2, ..., nbr_split, "test_dataset"
          `names<-`( c(rep("",nbr_split),"test_dataset") )
        
        return(out)
      }) %>%
      c(nui[-(1:nbr_split)])
  }
  
  
  # crossfit each method
  lapply(
    seq_along(Balancing_Method),
    function(met_idx) {
      met = Balancing_Method[[met_idx]]
      output = lapply(
        1:nbr_split,
        function(split_idx) {
          # determine the index of the splits used to train
          # the nuisance functions for the current method
          nui_split_idx = method_info[[met_idx]][[split_idx]][["nui_train_split_idx"]]
          
          # create a list containing the arguments needed to compute
          # the method
          input_nui = lapply(
            seq_along(nui_split_idx), 
            function(i) {
              method_info[[met_idx]][["nui_unique_nuisance_idx"]][i] %>%
                # get the nuisance function
                unique_nuisance[[.]] %>%
                # get the predicted values of the i-th nuisance function
                { .[[nui_split_idx[i]]][["predicted.value"]][["test_dataset"]] }
            }) %>%
            # rename according to the nuisance function name
            `names<-`(names(met[["nuisance_fct"]]) %>%
                        `[`( . %in% formalArgs(met[["method"]]) )
            )
          
          c(input_nui, test_dataset) %>%
            # removed unused args
            `[`( names(.) %in% formalArgs(met[["method"]]) ) %>%
            { try(do.call(met[["method"]], .)) }
        })
      if (any(sapply(output, inherits, "try-error"))) return(fail_output)
      
      output %>%
        do.call(cbind, .) %>%
        # compute the mean for each method
        apply(1, mean)
      
    }) %>%
    # set the names
    setNames(names(Balancing_Method))
}


crossfit = function(X_tr, Y_tr, Z_tr, X_te, Y_te, Z_te, Balancing_Method,
                    default_nbr_iter=100, default_nbr_fail=Inf) {
  
  # set the default value if not specified
  for (i in seq_along(Balancing_Method)) {
    # default for nbr_iter
    if (is.null(Balancing_Method[[i]]$nbr_iter)) {
      Balancing_Method[[i]]$nbr_iter = default_nbr_iter
    }
    
    # default for max_fail
    if (is.null(Balancing_Method[[i]]$max_fail)) {
      Balancing_Method[[i]]$max_fail = ifelse(
        is.null(default_nbr_fail),
        Balancing_Method[[i]]$nbr_iter,
        default_nbr_fail
      )
    }
    
    # default for nbr_split
    if (is.null(Balancing_Method[[i]]$nbr_split)) {
      Balancing_Method[[i]]$nuisance_fct %>% length
    }
  }
  
  # sort method by "nbr_split"
  Sorted_Balancing_Method = Balancing_Method %>%
    split(x = .,
          f = sapply(., `[[`, "nbr_split")
    )
  
  # cross-fit the methods nbr_iter time
  lapply(
    Sorted_Balancing_Method,
    function(Bal_Met) {
      browser()
      max_fail = sapply(Bal_Met, `[[`, "max_fail")
      nbr_crossfit = sapply(Bal_Met, `[[`, "nbr_iter")
      nbr_split = Bal_Met[[1]]$nbr_split
      
      Bal_Met_plus = add_method_info(Bal_Met)
      output = lapply(Bal_Met, function(met) NULL)
      idx_met = seq_along(Bal_Met)
      nbr_success = nbr_fail = setNames(rep(0, length(Bal_Met)), names(Bal_Met))
      
      while (length(idx_met) > 0) {
        # compute method on new split
        nbr_try = min(nbr_crossfit[idx_met] - nbr_success[idx_met])
        out = replicate(
          nbr_try,
          simple_crossfit(
            X_tr = X_tr, Y_tr = Y_tr, Z_tr = Z_tr,
            X_te = X_te, Y_te = Y_te, Z_te = Z_te,
            nbr_split = nbr_split,
            Balancing_Method = subset_Balancing_Method(Bal_Met_plus, idx_met)),
          simplify = FALSE
        )
        
        # save values in output
        for (i in seq_along(idx_met)) {
          j = idx_met[i]
          out_clean = lapply(out, `[[`, names(Bal_Met)[j]) %>%
            lapply(function(x) if (inherits(x, "try-error")) NULL else x) %>%
            do.call(rbind, .)
          output[[j]] %<>% rbind(out_clean)
          nbr_success[j] = nbr_success[j] + if (is.null(out_clean)) 0 else nrow(out_clean)
          nbr_fail[j] = nbr_fail[j] + nbr_try - nrow(out_clean)
        }
        
        idx_met = seq_along(Bal_Met)[(nbr_success < nbr_crossfit) & (nbr_fail < max_fail)]
      }
      
      # compute the median
      output[nbr_success == nbr_crossfit] %<>% lapply(function(met) apply(met, 2, median))
      output[nbr_success != nbr_crossfit] %<>% lapply(function(met) NA)
      
      # compute the median
      output = lapply(output, function(met) apply(as.matrix(met), 1, median))
      c(output, 'nbr_fail' = list(nbr_fail))
    }) %>%
    # concatenate and rename results
    {
      temp = do.call(c, .)
      names(temp) = lapply(., names) %>%
        do.call(c,.)
      temp
    }
}