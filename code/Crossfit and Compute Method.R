# This script implement the functions used to compute the crossfited estimator
# and non crossfited estimator.
# User should use the function "compute_method"

#### Initialization ####

debug_crossfit_and_comput = FALSE
if (debug_crossfit_and_comput) {
  setwd("G:/Mon Drive/Stage M2/Code")
  rm(list = ls()[ls() != "debug_crossfit_and_comput"])
  graphics.off()
  set.seed(123)
}
# debug
debug_var = 0

# Library import

library(magrittr) # pipeline



#### Functions ####


# get the index of the rows containing NA value
.idx_na_row = function(l) apply(
  l, 1, function(row) any(is.na(row)) ) %>%
  unlist %>% which

do.call_if_err_na = function(func, args) {
    try(do.call(func, args[names(args) %in% formalArgs(func)])) %>%
    # return NA if func fail
    { if (inherits(., "try-error")) NA else . }
}


check_method = function(Balancing_Method) {
  # check if all arguments / parameters are given to compute the method
  sapply(
    seq_along(Balancing_Method),
    function(met_idx)
    {
      met = Balancing_Method[[met_idx]]
      met_name = names(Balancing_Method)[met_idx]
      
      # names of the args which are going to be passed on to the function
      met_args_name = c("data") %>%
      { .[. %in% formalArgs(met$main)] } %>%
        c(names(met$aux))
      
      # check if their are enough args to run the function
      msg_met1 = formals(met$main) %>%
        # select args with no default value
      { names(.)[sapply(.,function(s) ifelse(is.symbol(s),s=="",FALSE))] } %>% 
      { . [ !(. %in%  met_args_name) ] } %>%
        paste0(collapse = '", "') %>%
        paste0('"', ., '"') %>%
        {ifelse(test = all(.=='""'),
                yes = "",
                no =paste0(
                  "\n", met_name, " : in \"main\", argument ", ., " is missing, with no default"
                ))}
      
      # check if all the args are used by the function
      msg_met2 = met_args_name %>%
      { .[ !(. %in% formalArgs(met$main)) ] } %>%
        paste0(collapse = '", "') %>%
        paste0('"', ., '"') %>%
        {
          ifelse(test = all(. == '""'),
                 yes = "",
                 no = paste0("\n", met_name, " : in \"aux\", unused element ", .)
          )
        }
      
      # same things as above but for the nuisance function(s) train
      msg_nui = sapply(
        met$aux,
        function(nui)
        {
          # check for the train function
          nui_train_args_name = c("data") %>%
          { .[. %in% formalArgs(nui$train)] }
          
          msg_train1 = formals(nui$train) %>%
            # select args without default value
          { names(.)[sapply(.,function(s) ifelse(is.symbol(s),s=="",FALSE))] } %>% 
          { .[ !(. %in% nui_train_args_name) ] } %>%
            paste0(collapse = '", "') %>%
            paste0('"', ., '"') %>%
            {
              ifelse(test = all(. == '""'),
                     yes = "",
                     no = paste0(
                       "\n ", names(met$aux)[met$aux %in% list(nui)],
                       " training function in ", met_name, " : argument ", .,
                       " is missing, with no default"
                     )
              )
            }
          
          msg_train2 = nui_train_args_name %>%
          { .[ !(. %in% formalArgs(nui$train)) ] } %>%
            paste0(collapse = '", "') %>%
            paste0('"', ., '"') %>%
            {
              ifelse(
                test = all(. == '""'),
                yes = "",
                no = paste0("\n ", names(met$aux)[met$aux %in% list(nui)],
                            " training function in ", met_name, " : unused argument ", .)
              )
            }
          
          # check for the predict function
          
          nui_pred_args_name = c("data", "trained.model") %>%
          { .[. %in% formalArgs(nui$predict)] }
          
          msg_pred1 = formals(nui$predict) %>%
            # select args without default value
          { names(.)[sapply(., function(s) ifelse(is.symbol(s),s=="",FALSE))] } %>%
          { .[ !(. %in% nui_pred_args_name) ] } %>%
            paste0(collapse = '", "') %>%
            paste0('"', ., '"') %>%
            {
              ifelse(test = all(. == '""'),
                     yes = "",
                     no = paste0(
                       "\n", met_name, " : in \"aux\" element \"",
                       names(met$aux)[met$aux %in% list(nui)], "\" predict function, argument \"",
                       ., "\" is missing, with no default"
                       # "\n", names(met$aux)[met$aux %in% list(nui)],
                       # " predict function in ", met_name, " : argument ", .,
                       # " is missing, with no default"
                     )
              )
            }
          
          msg_pred2 = nui_pred_args_name %>%
          { .[ !(. %in% formalArgs(nui$predict)) ] } %>%
            paste0(collapse = '", "') %>%
            paste0('"', ., '"') %>%
            {
              ifelse(test = all(. == '""'),
                     yes = "",
                     no = paste0(
                       "\n", names(met$aux)[met$aux %in% list(nui)],
                       " predict function in ", met_name, " : unused argument ", .)
              )
            }
          
          paste0(msg_train1, msg_train2, msg_pred1, msg_pred2) %>%
            return
        }) %>%
        paste0(collapse = "")
      
      
      msg = paste0(msg_met1, msg_met2, msg_nui)
      
      # check if the number of split is enough when specified
      
      if (is.null(met$crossfit)) return(msg)
      if (!is.logical(met$crossfit)) msg %<>% paste0("\n ", met_name, ' : "crossfit" must be logical or NULL.')
      if (!met$crossfit) return(msg)
      if (is.null(met$nbr_split)) return(msg)
      if (!is.numeric(met$nbr_split)) msg %<>% paste0("\n ", met_name, ' : "nbr_split" must be numeric or NULL')
      if ((length(met$aux) + 1) <= met$nbr_split) return(msg)
      
      paste0(msg, "\n", met_name, " :  Not enough split to crossfit the method") %>%
        return
    }) %>%
    {
      if (all(.=="")) "" else paste0(.[.!=""],collapse = "\n")
    } %>%
    cat()
}


non_crossfit = function(X, Y, Z, Balancing_Method, ...) {
  # Compute the Balancing_Method without cross fiting
  n = length(Z) # number of observation
  p = ncol(X) # number of covariate
  input = list("X" = X, "Y" = Y, "Z" = Z)
  
  # for each unique nuisance function, train the nuisance function and predict
  # the values on the same data set
  
  unique_nuisance = lapply(Balancing_Method, `[[`, "nuisance_fct") %>% 
    # concatenate the nuisance functions
    do.call(c,.) %>% 
    # remove duplicate
    unique %>%
    # train and predict the values on the data set
    lapply(
      function(nui)
      {
        input %>%
          # remove unused args
        { .[names(.) %in% formalArgs(nui$train)] } %>%
        { try(do.call(nui$train, .)) } %>%
          # return NA if "nui" fail
        { if (inherits(., "try-error")) NA else . } %>%
          c(input, "trained.model" = list(.)) %>%
          # remove unused args
          { .[names(.) %in% formalArgs(nui$predict)] } %>%
          # predict the values in the current split
          { try(do.call(nui$predict, .)) } %>%
          # return NA if "nui$predict" fail
          { if (inherits(., "try-error")) NA else . } %>%
          list %>%
          c(nui, "predicted.value" = .)
      })
  
  # compute the methods on each split
  
  lapply(
    Balancing_Method,
    function(met)
    {
      # create a list containing the args needed by "met$method"
      input_nui = lapply(
        seq_along(met$nuisance_fct), 
        function(i)
        {
          unique_nuisance %>%
            # select the train functions
            lapply(`[`, c("train","predict") ) %>%
            # get the index of the matching function
            match(met$nuisance_fct[i], .) %>%
            # get the preidcted values
            { unique_nuisance[[.]][["predicted.value"]] }
        })
      
      # rename elements of "input_nui"
      names(input_nui) = names(met$nuisance_fct)
      
      c(input, input_nui) %>%
        # remove unused args
      { .[names(.) %in% formalArgs(met$method)] } %>%
      { try(do.call(met$method, .)) } %>%
        #### test ####
        # return NA if "met$method" fail
      # { if (inherits(., "try-error")) NA else . } %>%
        unlist
    })
}


simple_crossfit = function(X, Y, Z, Balancing_Method, nbr_split) {
  n = length(Z) # number of observations
  p = ncol(X) # number of covariates
  
  # split the data evenly without overlap
  
  splited_data = cbind(X, Y, Z)[sample(n),] %>% # shuffle the rows
    # split [X, Y, Z] into "nbr_split" splits
    split( rep_len(1:nbr_split, n) ) %>% 
    # unbind [X, Y, Z]
    lapply( function(data) list(
      "X" = data[, 1:p],
      "Y" = data[, p+1],
      "Z" = data[, p+2])
    )
  
  
  # for each unique nuisance function, on each split train the nuisance
  # function and use the fitted model to predict the values on the
  # other splits
  
  unique_nuisance = lapply(Balancing_Method, `[[`, "nuisance_fct") %>% 
    # concatenate the nuisance functions
    do.call(c,.) %>% 
    # remove duplicate
    unique %>%
    lapply(
      function(nui) {
        lapply(
          1:nbr_split,
          function(i) {
            # train the nuisance function on split "i"
            
            trained.model = splited_data[[i]] %>%
              do.call_if_err_na(nui$train, .) %>%
              list()
            
            # use "trained.model" to predict the value on the other splits
            
            predicted.value = lapply(
              1:nbr_split,
              function(j) {
                # do not predict value on the same data set used to train the
                # nuisance function
                
                if (i == j) return(NA)
                
                # predict the values on split "j"
                
                c(splited_data[[j]], "trained.model" = trained.model) %>%
                  do.call_if_err_na(nui$predict, .)
              })
            
            list("trained.model"= trained.model,
                 "predicted.value" = predicted.value
            )
          }) %>%
          c(nui)
      })
  
  # crossfit each method
  
  lapply(
    Balancing_Method,
    function(met) {
      lapply(
        1:nbr_split,
        function(split_idx) {
          # determine the index of the splits used to train the
          # nuisance functions for the current split
          
          nui_split_idx = met$nuisance_fct %>%
            # get the minmum number of split to crossfit the method
            { length(.) + 1 } %>%
            # create an array centered around "split_idx"
            seq(from = split_idx - . %/% 2, length.out = .) %>%
            # remove the index of the current split
            { .[. != split_idx] } %>%
            # ensure the array is in [1, nbr_split]
            { (.-1) %% nbr_split + 1 }
          
          
          # create a list containing the arguments to use for the current method
          
          input_nui = lapply(
            seq_along(nui_split_idx), 
            function(i) {
              unique_nuisance %>%
                # select the train functions
                lapply(`[`, c("train","predict") ) %>%
                # get the index of the matching function
                match(met$nuisance_fct[i], .) %>%
                # get the predicted values
                { unique_nuisance[[.]][[nui_split_idx[i]]][["predicted.value"]][[split_idx]] }
            })
          
          # rename elements of "input_nui"
          names(input_nui) = names(met$nuisance_fct)
          
          c(input_nui, splited_data[[split_idx]]) %>%
            do.call_if_err_na(met$method, .) %>%
            unlist()
        }) %>%
        do.call(rbind, .) %>%
        # compute the mean for each method
        apply(2,mean)
    })
}


crossfit = function(X, Y, Z, Balancing_Method, default_nbr_crossfit=100, default_nbr_fail=NULL, ...) {
  # set the default value if not specified
  
  for (i in seq_along(Balancing_Method)) {
    # default for nbr_crossfit
    if (is.null(Balancing_Method[[i]]$nbr_crossfit)) {
      Balancing_Method[[i]]$nbr_crossfit = default_nbr_crossfit
    }
    
    # default for max_fail
    if (is.null(Balancing_Method[[i]]$max_fail)) {
      Balancing_Method[[i]]$max_fail = ifelse(
        is.null(default_nbr_fail),
        Balancing_Method[[i]]$nbr_crossfit,
        default_nbr_fail
      )
    }
    
    # default for nbr_split
    if (is.null(Balancing_Method$nbr_split)) {
      Balancing_Method[[i]]$nbr_split =
        Balancing_Method[[i]]$nuisance_fct %>% length + 1
    }
  }
  
  
  # sort method by "nbr_split"
  Sorted_Balancing_Method = Balancing_Method %>%
    split(x = .,
          f = sapply(., `[[`, "nbr_split"))
  
  
  # cross-fit the methods nbr_crossfit time
  lapply(
    Sorted_Balancing_Method,
    function(Bal_Met) {
      max_fail = sapply(Bal_Met, `[[`, "max_fail")
      nbr_crossfit = sapply(Bal_Met, `[[`, "nbr_crossfit")
      nbr_split = Bal_Met[[1]]$nbr_split
      
      output = lapply(Bal_Met, function(met) NULL)
      idx_met = seq_along(Bal_Met)
      nbr_success = nbr_fail = sapply(Bal_Met, function(met) 0)
      
      while (length(idx_met) > 0) {
        # compute method on new split
        nbr_try = min(nbr_crossfit[idx_met] - nbr_success[idx_met])
        out = replicate(
          nbr_try,
          simple_crossfit(
            X = X, Y = Y, Z = Z,
            nbr_split = nbr_split,
            Balancing_Method = Bal_Met[idx_met]
          ),
          simplify = FALSE
        )
        
        # save values in output
        for (i in seq_along(idx_met)) {
          j = idx_met[i]
          out_clean = lapply(out, `[[`, names(Bal_Met)[j]) %>%
            lapply(function(x) if (any(is.na(x))) NULL else x) %>%
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
      
      c(output, 'nbr_fail' = list(nbr_fail))
    }) %>%
    # concatenate and rename results
    {
      nbr_fail = lapply(., `[[`, 'nbr_fail') %>%
        {setNames(do.call(c, .), do.call(c, lapply(., names)))}
      temp = do.call(c, lapply(., function(met) met[-length(met)]))
      names(temp) = sapply(lapply(., function(met) met[-length(met)]), names) %>%
        do.call(c, .)
      c(temp, 'nbr_fail' = list(nbr_fail))
    }
}


compute_method = function(X, Y, Z, Balancing_Method, ...) {
  # compute the non crossfited method
  
  crossfit_bool = sapply(
    Balancing_Method,
    function(met) ifelse(is.null(met$crossfit), FALSE, met$crossfit)
  )
  
  non_crossfit_output = non_crossfit(
    X = X, Y = Y, Z = Z, Balancing_Method = Balancing_Method[!crossfit_bool], ...)

  if (!any(crossfit_bool)) return(non_crossfit_output)
  
  # compute the crossfited method
  crossfit_output =  crossfit(
    X = X, Y = Y, Z = Z, Balancing_Method = Balancing_Method[crossfit_bool], ...)
  
  c(non_crossfit_output, crossfit_output)
}


#### Test Functions ####

if (debug_crossfit_and_comput)
{
  # variables declaration
  
  met0 = function(X,Y,Z,nui_0){ list('ATE'=1,'ATT'=10) }
  met1 = function(X,Y,Z,nui_0,nui_2){ list('ATE'=-1, 'ATT'=-10) }
  nui0_train = function(X,Y,Z){return(0)}
  nui0_pred = function(X,Y,Z,trained.model){return(0)}
  nui1_train = function(X,Y,Z){return(1)}
  nui1_pred = function(X,Y,Z,trained.model){return(1)}
  
  
  method_a = list(
    
    'method' = met0,
    
    'nuisance_fct' = list(
      
      'nui_0' = list(
        'train' = nui0_train,
        'predict' = nui0_pred)
    ),
    
    'crossfit' = TRUE,
    'nbr_crossfit' = 10
  )
  
  method_c = list(
    
    'method' = met0,
    
    'nuisance_fct' = list(
      
      'nui_0' = list(
        'train' = nui0_train,
        'predict' = nui0_pred)
    ),
    
    'crossfit' = TRUE,
    'nbr_crossfit' = 3,
    'nbr_split' = 5
  )
  method_d = list(
    
    'method' = met1,
    
    'nuisance_fct' = list(
      
      'nui_0' = list(
        'train' = nui0_train,
        'predict' = nui0_pred),
      
      'nui_2' = list(
        'train' = nui1_train,
        'predict' = nui1_pred)
    ),
    
    'crossfit' = TRUE,
    'nbr_crossfit' = 10
  )
  
  method_b = list(
    
    'method' = met1,
    
    'nuisance_fct' = list(
      
      'nui_0' = list(
        'train' = nui0_train,
        'predict' = nui0_pred),
      
      'nui_2' = list(
        'train' = nui1_train,
        'predict' = nui1_pred)
    )
  )
  Balancing_Method = list("method A" = method_a,
                          "method B" = method_b,
                          "method C" = method_c,
                          "method D" = method_d)
  
  set.seed(123)
  
  n = 12
  X = data.frame(x1=1:n, x2=as.numeric(runif(n)<0.5))
  Y = runif(n) < 0.25
  Z = runif(n) < 0.75
  non_crossfit(X = X, Y = Y,Z = Z, Balancing_Method = Balancing_Method)
  crossfit(X, Y, Z, Balancing_Method,nbr_crossfit = 100)
  compute_method(X,Y,Z, Balancing_Method)
} else rm(debug_crossfit_and_comput)