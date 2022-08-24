#' R6 Class representing concentration-response data belonging to a compound
#'
#' Data belongig to a compound has a source file with a seperator, compound name, data and a protocol.
#' @export
CRDataTable <- R6::R6Class("CRDataTable",
                           private = list(
                             #' @description
                             #' Defining the average function used in this CRDataTable.
                             average_method = function(x){
                               
                               if(self$avg_method == 'mean'){
                                 return(as.double(mean(x, na.rm = T)))
                               }
                               
                               return(as.double(median(x, na.rm = T)))
                             },
                             
                             
                             #' @description
                             #' A generic helper function to apply filter on dt according m. Mainly used in methods to ensure data quality.                         
                             apply_filter = function(m, on_endpoints){
                               experiments <- m[, get('Experiment_ID')]
                               for(i in 1:length(experiments)){
                                 experiment_id <- experiments[i]
                                 for(j in 1:length(on_endpoints)){
                                   endpoint <- on_endpoints[j]
                                   if(m[get('Experiment_ID') == experiment_id, get(endpoint)]){
                                     self$dt[get('Experiment_ID') == experiment_id, (endpoint) := NA]
                                     paste0(experiment_id, endpoint, 'removed.') %push% self # log this process
                                   }
                                 }
                               }
                               'Filter Mask applied.'
                             },
                             
                             
                             #' @description
                             #' A generic helper function to add new column given by quotiens according counter_col and denominator_col. Mainly used in methods to generate endpoints induced by quotients.     
                             add_relative_col = function(counter_col, denominator_col, new_col = NA, delete_after=T){
                               if(! all(c(counter_col, denominator_col) %in% colnames(self$dt)) ){
                                 return('Cannot add relative column to data because specified counter or denomitor is not found in the endpoint columns.' %push% self)
                               }
                               
                               if(typeof(new_col) == 'character'){
                                 new_col_name_4_relative <- new_col
                               }
                               else{
                                 new_col_name_4_relative <- paste(counter_col, '/', denominator_col, sep = ' ')
                               }
                               
                               self$dt[, (new_col_name_4_relative) := 100*(get(counter_col)/get(denominator_col))]
                               # is.na(self$dt) <- sapply(self$dt, is.infinite)
                               self$dt[, (new_col_name_4_relative) := lapply(.SD, function(x) ifelse(is.infinite(x), 0, x)), .SD = new_col_name_4_relative] # x/0 case
                               self$dt[, (new_col_name_4_relative) := lapply(.SD, function(x) ifelse(is.nan(x), 0, x)), .SD = new_col_name_4_relative] # 0/0 case
                               
                               if(delete_after) {
                                 self$dt[,(counter_col):=NULL]
                               }
                               
                               'Add relative column sucessfully.' %push% self
                             },
                             
                             
                             #' @description
                             #' Normalization, Additive or Mulitplicative
                             #' @param endpt_conc_table data.table of with columns (Endpoint, Reference Concentration).
                             #' @param composition binary arithmetic operation, additive or multiplicative
                             normalize = function(endpt_conc_table, composition='multiplicative'){
                               N = nrow(endpt_conc_table)
                               endpoints <- self$get_endpoints()
                               for(i in 1:N){
                                 endpoint <- endpt_conc_table[i, get('Endpoint')]
                                 ref_conc <- endpt_conc_table[i, get('Concentration')]
                                 
                                 if(!(endpoint %in% endpoints)) next;
                                 
                                 # calc ctl averages for the endpoint and add helper col to self$dt
                                 averaged <- self$dt[get('Concentration') == ref_conc, lapply(.SD, private$average_method), .SDcols = endpoint, keyby = c('Experiment_ID')]
                                 tmp_ctl_col <- paste0('.Ctl.', endpoint)
                                 colnames(averaged)[colnames(averaged) %in% c(endpoint)] <- tmp_ctl_col
                                 self$dt <- self$dt[averaged, on = c('Experiment_ID')]
                                 
                                 # normalize endpoint
                                 if(composition == 'multiplicative'){
                                   self$dt[, (endpoint):= (100* get(endpoint)) / get(tmp_ctl_col)]
                                 }
                                 else if(composition == 'additive'){
                                   self$dt[, (endpoint):= get(endpoint) - get(tmp_ctl_col)]
                                   self$dt[, (endpoint) := lapply(.SD, function(x) ifelse(x<0, 0, x)), .SD = endpoint]
                                 }
                                 
                                 # delete helper
                                 self$dt[, c(tmp_ctl_col):= NULL]
                                 
                                 # log step
                                 paste(composition, 'normalization applied for', endpoint, 'relative to', ref_conc, '.', sep = ' ') %push% self
                               }
                               
                               # log step
                               'Normalization finished.' %push% self
                             },
                             
                             
                             #' @description
                             #' Function for Curve-Fitting that relies on the drc-package.
                             trydrm = function(d, e, respect_visual_criteria = TRUE){
                               result <- list('model' = NA, 'drcObject' = NA, 'aic' = NA)
                               general_model <- tryCatch({
                                 drm(endpoint ~ dose, data = d, robust = 'mean', fct = L.3())
                               },
                               error = function(cond){
                                 paste('error in initial fitting for', e, sep=' ') %push% self
                                 print(paste('error in initial fitting for', e, sep=' '))
                               })
                               
                               if(typeof(general_model) != 'list'){
                                 result$origData <- d
                                 return(result)
                               }
                               
                               # when initial fit suceeded, find accurate model and fit with that
                               models <- mselect(general_model, list(LL2.3(), LL2.4(), W1.3(), W1.4(), W2.3(), W2.4(), EXD.3(), EXD.2(), LL2.5(), BC.4(), BC.5()), linreg = F)
                               
                               flag <- TRUE
                               for(i in 1:length(rownames(models))){
                                 best_model <- rownames(models)[i]
                                 
                                 best_model <- paste0(best_model, '()')
                                 result$model <- best_model
                                 best_model <- eval(parse(text=best_model))
                                 
                                 # mselect does not ensure that the optimization performed in the drm function will terminate ! 
                                 result$drcObject <- try(drm(endpoint ~ dose, data = d, robust = 'mean', fct = best_model))
                                 if(typeof(result$drcObject) != 'list'){
                                   next
                                 }
                                 
                                 result$aic <- models[i, 'IC']
                                 result$origData <- result$drcObject$origData
                                 xlim_start = min(result$origData[["dose"]])
                                 xlim_end <- max(result$origData[["dose"]])*self$extrapolation
                                 
                                 neue_daten <- exp(seq(log(xlim_start), log(xlim_end), length = self$grid_size))
                                 neue_daten[1] <- xlim_start
                                 neue_daten[self$grid_size] <- xlim_end
                                 neue_daten <- data.frame(
                                   dose = neue_daten
                                 )
                                 
                                 p <- predict(result$drcObject, newdata =  neue_daten, interval = 'confidence')
                                 p <- cbind(neue_daten, p)
                                 colnames(p) <- c('dose', 'best', 'lower', 'upper')
                                 
                                 orig_doses <- result$origData[['dose']]
                                 orig_doses <- sort(unique(orig_doses))[3]
                                 pd <- as.data.table(p)
                                 pd <- na.omit(pd[, (c('dose')):= lapply(.SD, function(x) ifelse(x > orig_doses, x, NA)), .SDcols = c('dose')])
                                 pd <- as.data.frame(pd)
                                 
                                 
                                 # respect visual criteria to find a plausible fit
                                 if(respect_visual_criteria){
                                   if( (all(! is.na(p['lower']))) & (all(! is.na(p['upper']))) ){ # criterion 1
                                     if( all(p['lower'] < p['upper']) ) { # criterion 2
                                       if( all(pd['lower'] > 0) ){ # criterion 3
                                         flag <- FALSE
                                         break
                                       }
                                     }
                                   }
                                 }
                                 else{
                                   flag <- FALSE
                                   break
                                 }
                               }
                               
                               # if all possible models do not fulfill visual criteria
                               if(flag){
                                 for(i in 1:length(rownames(models))){
                                   best_model <- rownames(models)[i]
                                   
                                   best_model <- paste0(best_model, '()')
                                   result$model <- best_model
                                   best_model <- eval(parse(text=best_model))
                                   
                                   # mselect does not ensure that the optimization performed in the drm function will terminate ! 
                                   result$drcObject <- try(drm(endpoint ~ dose, data = d, robust = 'mean', fct = best_model))
                                   if(typeof(result$drcObject) != 'list'){
                                     next
                                   }
                                   
                                   result$aic <- models[i, 'IC']
                                   result$origData <- result$drcObject$origData
                                   
                                   break
                                 }
                                 
                                 result$model <- 'flagged'
                               }
                               
                               
                               # worst case - out of loop without any fit from mselect for some reason 
                               if(typeof(result$drcObject) != 'list'){
                                 paste('worst case - out of loop without any fit from mselect for', e, sep=' ') %push% self
                                 print(paste('worst case - out of loop without any fit from mselect for', e, sep=' '))
                                 result$drcObject <- general_model
                               }
                               
                               
                               
                               # linear-model vs result-fit
                               # if(FALSE){
                               lin_fit <- lm(data = d, formula = endpoint ~ dose)
                               lin_frame <- list('model' = 'lm', 'drcObject' = lin_fit, 'aic' = AIC(lin_fit), 'origData' = lin_fit[['model']])
                               N_dt = result$drcObject$sumList$lenData
                               
                               if(AIC(lin_fit) <= result$aic){
                                 if(!self$no_effect(d, e)){ # linear-model vs one-parameter-model
                                   result$model <- 'lm'
                                   result$drcObject <- lin_fit # actually not an drcObject
                                   result$aic <- AIC(lin_fit)
                                   result$origData <- lin_fit[['model']]
                                 }
                                 else{
                                   one_param_fit <- lm(data = d, formula = endpoint ~ 1)
                                   result$model <- 'no-effect'
                                   result$drcObject <- one_param_fit # actually not an drcObject
                                   result$aic <- AIC(one_param_fit)
                                   result$origData <- d
                                 }
                               }
                               # }
                               
                               
                               
                               result
                             }
                             
                           ),
                           public = list(
                             #' @field source_file is the path to source *.csv-file.
                             source_file = NULL,
                             
                             #' @field output folder name.
                             output_folder = NULL,
                             
                             #' @field sep seperator for *.csv-file.
                             sep = NULL,
                             
                             #' @field avg_method
                             avg_method = NULL,
                             
                             #' @field compound name of compound.
                             compound = NULL,
                             
                             #' @field dt raw data in form of a data.table.
                             dt = NULL,
                             
                             #' @field key_process counter.
                             key_process = NULL,
                             
                             #' @field solvent_control wording
                             solvent_control = NULL,
                             
                             #' @field renormalized table
                             renormalized = NULL,
                             
                             #' @field grid_size to use for curves
                             grid_size = NULL,
                             
                             #' @field extrapolation factor (for concentration range)
                             extrapolation = NULL,
                             
                             #' @field bg_normalized table 
                             bg_normalized = NULL,
                             
                             #' @field ctl_normalized marker
                             ctl_normalized = NULL,
                             
                             #' @field num_trt_ok requirements treatment groups
                             num_trt_ok = NULL,
                             
                             #' @field num_ctl_rpl_ok requirements control replicates
                             num_ctl_rpl_ok = NULL,
                             
                             #' @field cut_off_appl boolean whether custom cut off is applied
                             cut_off_appl = NULL,
                             
                             #' @field min_num_trt number of treatment groups
                             min_num_trt = NULL,
                             
                             #' @field min_num_ctl_rpl number of ctl repl
                             min_num_ctl_rpl = NULL,
                             
                             #' @field yaxis_lim limit for plots
                             yaxis_lim = NULL,
                             
                             #' @field protocol vector of characters.
                             protocol = c(),
                             
                             #' @field crc_dt boolean 
                             crc_dt = NULL,
                             
                             #' @field bmcs data.frame of the bmcs.
                             bmcs = data.frame(),
                             
                             #' @description
                             #' Create a new data belongig to a compound object.
                             #' @param source_file path to source *.csv-file.
                             #' @param ddata data.table of data, when data of source_file is available direct in that form 
                             #' @param sep seperator for *.csv-file.
                             #' @param compound compound name
                             #' @param general_setting data.frame with general settings
                             #' @return A new `CRDataTable` object.
                             initialize = function(source_file = NA, ddata, sep = ',', compound = NA, general_setting = NA){
                               'Create object' %push% self
                               self$output_folder <- 'output'
                               self$key_process <- 0
                               self$source_file <- source_file
                               self$compound <- compound
                               self$sep <- sep
                               self$bmcs <- data.frame(endpoint = character(), bmr = double(), lower_bmc = double(), bmc = character(), upper_bmc = character(), lbmc_boot = double(), bmc_boot = double(), ubmc_boot = double())
                               self$renormalized <- data.table('Endpoint' = character())
                               self$bg_normalized <- FALSE
                               self$ctl_normalized <- FALSE
                               self$num_trt_ok <- FALSE
                               self$num_ctl_rpl_ok <- FALSE
                               self$cut_off_appl <- FALSE
                               
                               if(is.na(source_file)){
                                 self$dt <- ddata
                               }
                               else{
                                 self$dt <- tryCatch({
                                   tmp <- fread(source_file, sep = sep)
                                   
                                   'Successfully reading file.' %push% self
                                   
                                   if(all(c('Experiment_ID', 'Concentration', 'Well') %in% colnames(tmp))){
                                     tmp
                                   }
                                   else {
                                     'Requiered columns not contained in the file! check the seperator used in the file.' %push% self
                                     data.table()
                                   }
                                 }, 
                                 error=function(cond) {
                                   'Error reading file.' %push% self
                                   return(data.table())
                                 },
                                 warning=function(cond) {
                                   cond$message %push% self 
                                   return(data.table())
                                 }
                                 )                                 
                               }

                               self$crc_dt <- list()
                               
                               # load general settings
                               if(typeof(general_setting) == 'list'){
                                 self$solvent_control <- general_setting[get('Setting')=='Your wording for Solvent Control', get('Value')]
                                 self$avg_method <- general_setting[get('Setting')=='Averaging method', get('Value')]
                                 self$grid_size <- as.numeric(general_setting[get('Setting')=='Grid size', get('Value')])
                                 self$extrapolation <- as.numeric(general_setting[get('Setting')=='Extrapolation factor', get('Value')])
                                 self$min_num_trt <- as.numeric(general_setting[get('Setting')=='Minimal Number of Treatment Groups', get('Value')])
                                 self$min_num_ctl_rpl <- as.numeric(general_setting[get('Setting')=='Minimal Number of Control Replicates', get('Value')])
                                 self$yaxis_lim <- as.numeric(general_setting[get('Setting')=='Y-Axis Limit', get('Value')])
                                 
                                 # write.table(general_setting, file = paste0(self$output_folder, '/Settings.csv'), sep=';', row.names = FALSE)
                               }
                               else{
                                 self$solvent_control <- 'Solvent control (SC)'
                                 self$avg_method <- 'median'
                                 self$grid_size <- 10000
                                 self$extrapolation <- 1.5
                                 self$min_num_trt <- 5
                                 self$min_num_ctl_rpl <- 3
                                 self$yaxis_lim <- 200
                               }
                               
                               # convert entries to numeric
                               endpoints <- self$get_endpoints()
                               self$dt[, (endpoints):= lapply(.SD, as.numeric), .SDcols = endpoints]
                               
                               self$save_dt('Input')
                               'Data Object created.' %push% self
                             },
                             
                             
                             #' @description
                             #' Get concentrations.
                             #' @return List of concentrations in the `CRDataTable` object.                         
                             get_doses = function(){
                               return(sort(unique(self$dt[, get('Concentration')])))
                             },
                             
                             
                             #' @description
                             #' Get endpoints.
                             #' @return List of endpoints in the `CRDataTable` object.
                             get_endpoints = function(){
                               return(colnames(self$dt[, !c('Experiment_ID', 'Concentration', 'Well')]))
                             },
                             
                             
                             #' @description
                             #' Get experiments.
                             #' @return List of experiments in the `CRDataTable` object.                       
                             get_experiments = function(){
                               unique(self$dt[, get('Experiment_ID')])
                             },
                             
                             
                             #' @description
                             #' Add new endpoints by quotiens.
                             add_quotient_parameter = function(quotients = NA){
                               if(typeof(quotients) == 'list'){
                                 N = nrow(quotients)
                                 for(i in 1:N){
                                   numerator <- quotients[i, get('Numerator parameter')]
                                   denumerator <- quotients[i, get('Denominator parameter')]
                                   new_parameter <- quotients[i, get('Resulting parameter')]
                                   private$add_relative_col(counter_col = numerator, denominator_col = denumerator, new_col = new_parameter, delete_after = FALSE)
                                 }
                                 'New parameter given by quotients are added.' %push% self
                               }
                             },
                             
                             
                             #' @description
                             #' Export current state of dt as csv-file.
                             save_dt = function(foldername){
                               dir.create(paste0(self$output_folder, '/', self$key_process, '_', foldername))
                               write.table(self$dt, file = paste0(self$output_folder, '/', self$key_process, '_', foldername, '/', self$compound, '.csv'), sep=';', row.names = FALSE)
                               self$key_process <- self$key_process +1
                             },
                             
                             
                             #' @description
                             #' Apply custom cut off on dt according to cutoff_table.
                             custom_cut_off = function(cutoff_table = NA){
                               if(typeof(cutoff_table) == 'list'){
                                 N = nrow(cutoff_table)
                                 if(N > 0){
                                   endpoints <- self$get_endpoints()
                                   for(i in 1:N){
                                     coi <- cutoff_table[i, get('Concentration')]
                                     poi <- cutoff_table[i, get('Parameter of interest')]
                                     
                                     if(!(poi %in% endpoints)) next;
                                     
                                     cutoff <- as.numeric(cutoff_table[i, get('Cut off')])
                                     aff_par <- cutoff_table[i, get('Affected parameter')]
                                     
                                     averaged <- self$dt[get('Concentration') == coi, lapply( .SD, private$average_method), .SDcols = c(poi), keyby = c('Experiment_ID')]
                                     experiments <- averaged[, get('Experiment_ID')]
                                     
                                     for(j in 1:length(experiments)){
                                       experiment_id <- experiments[j]
                                       avg <- averaged[get('Experiment_ID') == experiment_id, get(poi)]
                                       if( is.na(avg) || avg <= cutoff )
                                         self$dt[get('Experiment_ID') == experiment_id, (aff_par) := NA]
                                     }
                                   }
                                   self$cut_off_appl <- TRUE
                                   self$save_dt('CustomCutoff')                                   
                                 }
                               }
                               'Custom cut-offs applied.' %push% self
                             },
                             
                             
                             #' @description
                             #' Ensure data quality according number of conditions.                                                  
                             conditions = function(){
                               endpoints <- self$get_endpoints()
                               
                               mask <- self$dt[, lapply(.SD, median, na.rm = T), .SDcols = endpoints, keyby = c('Experiment_ID', 'Concentration')]
                               mask <- mask[, lapply(.SD, function(v) length(na.omit(v)) < self$min_num_trt), .SDcols = endpoints, keyby = c('Experiment_ID')]
                               
                               private$apply_filter(mask, on_endpoints = endpoints)
                               
                               self$num_trt_ok <- TRUE
                               self$save_dt('RemovedExperimentsWithTooLessConditions')
                               'Removed experiments with too less conditions.' %push% self
                             },
                             
                             
                             #' @description
                             #' Ensure data quality according number of control (technical) replicates.
                             control_replicates = function(){
                               endpoints <- self$get_endpoints()
                               mask <- self$dt[get('Concentration')==self$solvent_control, lapply(.SD, function(v) length(na.omit(v)) < self$min_num_ctl_rpl), .SDcols = endpoints, keyby = c('Experiment_ID')]
                               private$apply_filter(mask, on_endpoints = endpoints)
                               
                               self$num_ctl_rpl_ok <- TRUE
                               self$save_dt('RemovedExperimentsWithTooLessControlReplicates')
                               'Removed experiments with too less control replicates.' %push% self
                             },
                             
                             
                             #' @description
                             #' Apply background normalization according to the background_table.                      
                             background_subtraction = function(background_table = NA){
                               if(typeof(background_table) == 'list'){
                                 private$normalize(background_table, composition = 'additive')
                                 
                                 self$bg_normalized <- TRUE
                                 self$save_dt('BackgroundNormalization')
                                 'Background Subtraction applied.' %push% self
                               }
                               else{
                                 'Background Subtraction not applied.' %push% self
                               }
                             },
                             
                             
                             #' @description
                             #' Apply control normalizaiton with respect to dr_endpoints_table.                                               
                             control_normalization = function(dr_endpoints_table = NA){
                               endpoints <- self$get_endpoints()
                               if(typeof(dr_endpoints_table) == 'list'){ # dynamic range normalization
                                 N = nrow(dr_endpoints_table)
                                 if(N > 0){
                                   for(i in 1:N){
                                     dr_endpoint <- dr_endpoints_table[i, get('Endpoint')]
                                     
                                     # skip if dr_endpoint is not avaible
                                     if(!(dr_endpoint %in% endpoints)) next;
                                     
                                     min_ctl <- dr_endpoints_table[i, get('Min control')]
                                     max_ctl <- dr_endpoints_table[i, get('Max control')]
                                     
                                     # calc min_ctl averages for the dr_endpoint and add helper col to self$dt
                                     averaged <- self$dt[get('Concentration') == min_ctl, lapply( .SD, private$average_method), .SDcols = c(dr_endpoint), keyby = c('Experiment_ID')]
                                     tmp_min_ctl_col <- paste0('.MinCtl.', dr_endpoint)
                                     colnames(averaged)[colnames(averaged) %in% c(dr_endpoint)] <- tmp_min_ctl_col
                                     self$dt <- self$dt[averaged, on = c('Experiment_ID')]
                                     
                                     # calc max_ctl averages for the dr_endpoint and add helper col to self$dt
                                     averaged <- self$dt[get('Concentration') == max_ctl, lapply( .SD, private$average_method), .SDcols = c(dr_endpoint), keyby = c('Experiment_ID')]
                                     tmp_max_ctl_col <- paste0('.MaxCtl.', dr_endpoint)
                                     colnames(averaged)[colnames(averaged) %in% c(dr_endpoint)] <- tmp_max_ctl_col
                                     self$dt <- self$dt[averaged, on = c('Experiment_ID')]
                                     
                                     # normalize dr_endpoint
                                     self$dt[, (dr_endpoint):= (100*(get(tmp_max_ctl_col) - get(dr_endpoint))) / (get(tmp_max_ctl_col) - get(tmp_min_ctl_col))]
                                     self$dt[, (dr_endpoint) := lapply(.SD, function(x) ifelse(x<0, 0, x)), .SD = dr_endpoint]
                                     
                                     # delete helper
                                     self$dt[, c(tmp_min_ctl_col, tmp_max_ctl_col):= NULL]
                                     
                                     # log step
                                     paste('Dynamic Range Normalization applied for', dr_endpoint, '.', sep = ' ') %push% self
                                   }
                                   
                                   dr_endpoints <- dr_endpoints_table[, get('Endpoint')]
                                   dr_endpoints <- dr_endpoints[ (dr_endpoints %in% endpoints) ]
                                   
                                   endpoints <- endpoints[! (endpoints %in% dr_endpoints)]
                                   
                                   # log step
                                   'Dynamic Range Normalization finished.' %push% self 
                                 }
                               }
                               
                               endpt_conc_table <- data.table('Endpoint' = endpoints, 'Concentration' = self$solvent_control)
                               private$normalize(endpt_conc_table)
                               
                               # log step
                               self$ctl_normalized <- TRUE
                               self$save_dt('ControlNormalization')
                               'Control Normalization finished.' %push% self
                             },
                             
                             
                             #' @description
                             #' Apply renormalization according to drc models drc_models.
                             renormalize = function(experiments = c('merged', 'single'), endpoint, drc_models){
                               
                               if(endpoint %in% self$renormalized[, get('Endpoint')]){ # skip if already renormalized!
                                 return(paste('Renormalization of', endpoint, 'according to the fit-model already applied.', sep = ' ') %push% self)
                               }
                               
                               # parameter
                               experiments <- match.arg(experiments)
                               
                               all_experiments <- names(drc_models)
                               lapply(all_experiments, function(experiment_id){
                                 best_fit <- drc_models[[experiment_id]]
                                 
                                 if(is.na(best_fit$model)){ # no-fit case
                                   return()
                                 }
                                 
                                 drc_model <- best_fit$drcObject
                                 if(best_fit$model %in% c('lm', 'no-effect')){
                                   lowest_dose_predres <- drc_model[["fitted.values"]][[1]]
                                 }
                                 else{
                                   lowest_dose_predres <- drc_model[["predres"]][, 'Predicted values'][[1]]
                                 }
                                 
                                 if(experiments == 'single'){
                                   self$dt[get('Experiment_ID') == experiment_id, (endpoint) := (100*get(endpoint)) / lowest_dose_predres ]
                                 }
                                 else if(experiments == 'merged'){ # in this case all_experiments = ['merged'] 
                                   self$dt[, (endpoint) := (100*get(endpoint)) / lowest_dose_predres ]
                                 }
                               })
                               
                               # update status of self$dt
                               tmp <- data.table('Endpoint' = endpoint)
                               self$renormalized <- rbind(self$renormalized, tmp)
                               
                               # self$save_dt('Renormalized')
                               paste('Renormalization of', endpoint, 'according to the fit-model applied.', sep = ' ') %push% self
                             },
                             
                             
                             #' @description
                             #' Get data / model for endpoint with seperation of data regarding all plate (merged) or single plates (single)
                             #' @param sth select data or model.
                             #' @param experiments merged or single option.
                             #' @param for_endpoint character or list of endpoints. If not restricted, for all endpoints will be printed.
                             #' @return Named list of the models (list('model' = ..., 'drcObject' = ...)) of the best fit. The names are the plate ids (IUF: Experiment_IDs) for experiments option single, and 'merged' for experiments option merged.                         
                             to_get = function(sth = c('data', 'model'), experiments = c('merged', 'single'), for_endpoint){
                               print(paste0("**", for_endpoint, "**"))
                               experiments <- match.arg(experiments)
                               sth <- match.arg(sth)
                               
                               # prepare data for fitting via drc package
                               endpoints <- self$get_endpoints()
                               
                               # prepare dt_endpoint table for both cases
                               if(experiments == 'merged'){
                                 averaged <- self$dt[, lapply(.SD, private$average_method), .SDcols = endpoints, keyby = c('Experiment_ID', 'Concentration')]
                                 cols <- c('Concentration', for_endpoint)
                                 dt_endpoint <- averaged[, ..cols]
                                 experiment_id_col <- data.table('Experiment_ID' = 'merged')
                                 dt_endpoint <- cbind(experiment_id_col, dt_endpoint)
                               }
                               else if(experiments == 'single'){
                                 cols <- c('Experiment_ID', 'Concentration', for_endpoint)
                                 dt_endpoint <- self$dt[, ..cols]
                               }
                               
                               # remove rows with non-numeric entries in column 'Concentration'
                               dt_endpoint[, ('Concentration') := as.numeric(get('Concentration'))]
                               dt_endpoint <- dt_endpoint[!is.na(get('Concentration'))] # na.omit
                               
                               all_experiments <- unique(dt_endpoint[, get('Experiment_ID')])
                               res <- lapply(all_experiments, function(experiment_id){
                                 cols <- c('Concentration', for_endpoint)
                                 dt_endpoint_experiment <- dt_endpoint[get('Experiment_ID') == experiment_id, ..cols]
                                 
                                 colnames(dt_endpoint_experiment)[1] <- 'dose'
                                 colnames(dt_endpoint_experiment)[2] <- 'endpoint'
                                 
                                 # if(length(na.omit(dt_endpoint_experiment[, endpoint])) == 0){
                                 # print("NO data")
                                 # return(NA)
                                 # }
                                 
                                 if(sth == 'model'){
                                   result <- private$trydrm(dt_endpoint_experiment, for_endpoint) # list('model'= , 'drcObject'=)
                                   
                                   if(experiment_id == 'merged'){
                                     averaged <- self$dt[, lapply(.SD, private$average_method), .SDcols = endpoints, keyby = c('Experiment_ID', 'Concentration')]
                                     cols_tmp <- c('Experiment_ID', 'Concentration', for_endpoint)
                                     averaged <- averaged[, ..cols_tmp]
                                     
                                     # remove rows with non-numeric entries in column 'Concentration'
                                     averaged[, ('Concentration') := as.numeric(get('Concentration'))]
                                     averaged <- averaged[!is.na(get('Concentration'))] # na.omit
                                     
                                     colnames(averaged)[1] <- 'experiment'
                                     colnames(averaged)[2] <- 'dose'
                                     colnames(averaged)[3] <- 'endpoint'
                                     
                                     result$origData <- averaged
                                   }
                                   else{
                                     tmp_dt <- self$dt
                                     cols_tmp <- c('Experiment_ID', 'Concentration', for_endpoint)
                                     tmp_dt <- tmp_dt[get('Experiment_ID') == experiment_id, ..cols_tmp]
                                     
                                     # remove rows with non-numeric entries in column 'Concentration'
                                     tmp_dt[, ('Concentration') := as.numeric(get('Concentration'))]
                                     tmp_dt <- tmp_dt[!is.na(get('Concentration'))] # na.omit
                                     
                                     colnames(tmp_dt)[1] <- 'experiment'
                                     colnames(tmp_dt)[2] <- 'dose'
                                     colnames(tmp_dt)[3] <- 'endpoint'
                                     
                                     result$origData <- tmp_dt
                                   }
                                   
                                   result
                                 }
                                 else if(sth == 'data'){
                                   dt_endpoint_experiment
                                 }
                               })
                               names(res) <- all_experiments
                               res
                             },
                             
                             
                             #' @description
                             #' Check if there is no-effect for given endpoint, data regarding given doses.
                             #' Therefore the nested models linear-regression and one-paramter-model are used, and the hypothesis is tested via tha likelihood-ratio test
                             #' implicitly using Wilks' theorem about asymptotic distribution of the ratio.  
                             #' @param experiments defines the seperation method.
                             #' @param for_endpint endpoint.
                             #' @param regarding list of doses.
                             #' @return Named list of boolean with names experiments.                      
                             no_effect = function(d, e, regarding = NA){
                               if(all(is.na(regarding))){
                                 regarding <- na.omit(d[, dose])
                               }
                               d <- d[dose %in% regarding]
                               d[d=='Inf'] <- NA
                               d <- na.omit(d)
                               if(all(is.na(d))){
                                 return(FALSE)
                               }
                               one_parameter_fit <- lm(data = d, formula = endpoint ~ 1)
                               lin_fit <- lm(data = d, formula = endpoint ~ dose)
                               
                               log_lik_null <- logLik(one_parameter_fit) # null hypothesis
                               log_lik_lf <- logLik(lin_fit)
                               lik_ratio <- -2*(log_lik_null - log_lik_lf)
                               deg_freedom_diff <- attr(log_lik_lf, "df") - attr(log_lik_null, "df")
                               p_value <- 1 - pchisq(lik_ratio, deg_freedom_diff)
                               return(!(p_value <= 0.05))
                             },
                             
                             
                             #' @description
                             #' Fit curves to the concentration response data.
                             #' @param experiments merged or single option.
                             #' @param for_endpoint character or list of endpoints. If not restricted, for all endpoints will be printed.
                             #' @param normal TRUE (default) or FALSE for renormalisation of the data regarding the smallest dose.
                             #' @return Plot of the data regarding to the endpoints (and experiments) with the benchmark-response and benchmark-dose.                         
                             fit_dt = function(experiments = c('merged', 'single'), endpoints_table = NA, rnrm = 'no'){
                               # parameter
                               experiments <- match.arg(experiments)
                               
                               if((typeof(endpoints_table) != 'list') || (nrow(endpoints_table) == 0) ){
                                 endpoints_table <- data.table('Name of PDF file' = gsub("/", "", self$get_endpoints()), 'Name in PDF'=self$get_endpoints(), 'Input'=self$get_endpoints(), 'Renormalization'=rnrm)
                               }
                               
                               N <- nrow(endpoints_table)
                               
                               for(i in 1:N) {
                                 endpoint <- endpoints_table[i, get('Input')]
                                 renormal <- ifelse(endpoints_table[i, get('Renormalization')] == 'yes', TRUE, FALSE) 
                                 if(!(endpoint %in% self$get_endpoints())) next;
                                 
                                 if(renormal) { # renormalization
                                   self$crc_dt[[endpoint]] <- self$to_get(sth = 'model', experiments = 'single', for_endpoint = endpoint)
                                   self$renormalize(experiments = 'single', endpoint, self$crc_dt[[endpoint]])
                                 }
                                 self$crc_dt[[endpoint]] <- self$to_get(sth = 'model', experiments = experiments, for_endpoint = endpoint)
                               }
                               
                               # log status
                               if(nrow(self$renormalized) > 0) self$save_dt('Renormalized');
                               
                               write.table(endpoints_table, file = paste0(self$output_folder, '/EndpointsReference.csv'), sep=';', row.names = FALSE)
                               'Concentration response data are fitted.' %push% self
                             },                    
                             
                             
                             #' @description
                             #' BMC calculation based on confidence bands (confidence interval projection) of the drc Package (in fact it is the delta method)
                             #' @param method Currently just the ci_pr option.
                             #' @param best_fit fit object.
                             #' @param bmr double or list of double with the benchmark-response
                             #' @return Plot of the data regarding to the endpoints (and experiments) with the benchmark-response and benchmark-dose.                         
                             bmc_calc = function(method = 'ci_pr', plate, for_endpoint, bmr){
                               
                               result <- c()
                               
                               
                               # bmc (estimation)
                               for(j in 1:length(bmr)){
                                 bmc_datarow <- c()
                                 
                                 bmr_j <- bmr[[j]]
                                 
                                 bmc_datarow$endpoint <- for_endpoint
                                 bmc_datarow$bmr <- 100 - bmr_j
                                 
                                 # visual way of estimation                                
                                 # lower bound of the confidence interval of the bmc estimation 
                                 translation_lower <- plate$lower - (100 - bmr_j)
                                 lower_bmc <- NA
                                 
                                 # bmc of the fit 
                                 translation <- plate$best - (100 - bmr_j)
                                 bmc <- NA
                                 
                                 # upper bound of the confidence interval of the bmc estimation 
                                 translation_upper <- plate$upper - (100 - bmr_j)
                                 upper_bmc <- NA
                                 
                                 translated <- data.table(dose = plate$dose, lower = translation_lower, best = translation, upper = translation_upper)
                                 translated$row_num <- seq.int(nrow(translated))
                                 
                                 # bmc candidates
                                 if(nrow(translated[lower <= 0]) > 0) {
                                   if(nrow(translated[lower >= 0]) > 0){
                                     
                                     
                                     lq <- translated[lower <= 0]
                                     gq <- translated[lower >= 0]
                                     
                                     lq <- lq[order(row_num)]
                                     for(i in 1:nrow(lq)){
                                       if((i+1) %in% gq[, get('row_num')]){
                                         lower_bmc <- lq[i, get('dose')]
                                         break
                                       }
                                     }
                                     
                                     
                                   }
                                 }
                                 
                                 if(nrow(translated[best <= 0]) > 0) {
                                   if(nrow(translated[best >= 0] > 0)){
                                     
                                     lq <- translated[best <= 0]
                                     gq <- translated[best >= 0]
                                     
                                     lq <- lq[order(row_num)]
                                     for(i in 1:nrow(lq)){
                                       if((i+1) %in% gq[, get('row_num')]){
                                         bmc <- lq[i, get('dose')]
                                         break
                                       }
                                     }
                                     
                                     
                                   }
                                 }
                                 
                                 if(nrow(translated[upper <= 0]) > 0) {
                                   if(nrow(translated[upper >= 0] > 0)){
                                     
                                     
                                     lq <- translated[upper <= 0]
                                     gq <- translated[upper >= 0]
                                     
                                     lq <- lq[order(row_num)]
                                     for(i in 1:nrow(lq)){
                                       if((i+1) %in% gq[, get('row_num')]){
                                         upper_bmc <- lq[i, get('dose')]
                                         break
                                       }
                                     }
                                     
                                   }
                                 }
                                 
                                 
                                 if(is.na(bmc)){
                                   lower_bmc <- NA
                                   upper_bmc <- NA
                                 }
                                 else{  # induction case
                                   if((!is.na(lower_bmc)) & (!is.na(upper_bmc))){
                                     if(lower_bmc > upper_bmc){
                                       abstellplatz <- upper_bmc
                                       upper_bmc <- lower_bmc
                                       lower_bmc <- abstellplatz
                                     }
                                   }
                                   else if(!is.na(lower_bmc)){
                                     if(lower_bmc > bmc){
                                       abstellplatz <- lower_bmc
                                       lower_bmc <- NA
                                       upper_bmc <- abstellplatz
                                     }
                                   }
                                   else if(!is.na(upper_bmc)){
                                     if(upper_bmc < bmc){
                                       abstellplatz <- upper_bmc
                                       upper_bmc <- NA
                                       lower_bmc <- abstellplatz
                                     }
                                   }
                                 }
                                 
                                 # save results
                                 bmc_datarow$lower_bmc <- lower_bmc
                                 bmc_datarow$bmc <- bmc
                                 bmc_datarow$upper_bmc <- upper_bmc
                                 
                                 self$bmcs <- rbind(self$bmcs, bmc_datarow)
                                 result <- rbind(result, bmc_datarow)
                               }
                               
                               result
                             },
                             
                             
                             #' @description
                             #' Print curves of the endpoints.
                             #' @param experiments merged or single option.
                             #' @param for_endpoint character or list of endpoints. If not restricted, for all endpoints will be printed.
                             #' @param normal TRUE (default) or FALSE for renormalisation of the data regarding the smallest dose.
                             #' @param bmr double or list of double with the benchmark-response
                             #' @return Plot of the data regarding to the endpoints (and experiments) with the benchmark-response and benchmark-dose.                         
                             plotc = function(experiments = c('merged', 'single'), for_endpoint = NULL, bmr = NA){
                               # parameter
                               experiments <- match.arg(experiments)
                               
                               if((typeof(for_endpoint) == 'list') && (nrow(for_endpoint) > 0)){ # no slashes allowed in pdfname column
                                 endpoints_final_names_table <- for_endpoint
                               }
                               else{
                                 endpoints_final_names_table <- data.table('Name of PDF file' = gsub("/", "", self$get_endpoints()), 'Name in PDF'=self$get_endpoints(), 'Input'=self$get_endpoints())
                               }
                               
                               
                               N <- nrow(endpoints_final_names_table)
                               
                               # bmr param                          
                               if(typeof(bmr) == 'list'){
                                 bmr <- as.list(bmr[, get('BMR')])
                               }
                               else {
                                 bmr <- list(0.0)
                               }
                               
                               result <- data.frame('dose' = double(), 'best' = double(), 'lower' = double(), 'upper' = double(), 'endpoint' = character())
                               
                               
                               
                               path <- paste0(self$output_folder, '/concentration-response-curves')
                               dir.create(path)
                               
                               path <- paste0(path, '/', self$compound)
                               dir.create(path)
                               
                               # for future
                               if(FALSE){
                                 curve_fits <- self$crc_dt
                                 if(length(curve_fits) == 0){
                                   self$fit_dt
                                 }
                               }
                               
                               # assumption correct/consistent mode (variable 'experiment')
                               curve_fits <- self$crc_dt
                               
                               
                               new_endpoints_names <- endpoints_final_names_table[, get('Name in PDF')]
                               
                               
                               group_colors <- c(rep('#00167B', length(new_endpoints_names)), rep('#B9B7BD', length(bmr)))
                               group_shapes <- c(rep(18, length(new_endpoints_names)), rep(NA, length(bmr)))
                               group_linetype <- c(rep('solid', length(new_endpoints_names)), rep('dashed', length(bmr)))
                               group_alpha <- c(rep(1, length(new_endpoints_names)), rep(.2, length(bmr)))
                               group_size <- c(rep(1.6, length(new_endpoints_names)), rep(.5, length(bmr)))
                               
                               
                               names(group_colors) <- c(new_endpoints_names, paste0('BMR', unlist(bmr)))
                               names(group_shapes) <- c(new_endpoints_names, paste0('BMR', unlist(bmr)))
                               names(group_linetype) <- c(new_endpoints_names, paste0('BMR', unlist(bmr)))
                               names(group_alpha) <- c(new_endpoints_names, paste0('BMR', unlist(bmr)))
                               names(group_size) <- c(new_endpoints_names, paste0('BMR', unlist(bmr)))
                               
                               
                               for(i in 1:N){
                                 endpoint <- endpoints_final_names_table[i, get('Input')]
                                 
                                 if(!(endpoint %in% self$get_endpoints())) next;
                                 
                                 new_endpoint_name <- endpoints_final_names_table[i, get('Name in PDF')]
                                 new_endpoint_pdf_name <- endpoints_final_names_table[i, get('Name of PDF file')]
                                 new_endpoint_pdf_name <- gsub("/", " ", new_endpoint_pdf_name) # to ensure valid file name of the output pdf
                                 # print(new_endpoint_name)
                                 
                                 grDevices::pdf(file = file.path(path, paste0(self$compound, new_endpoint_pdf_name, '.pdf')), 9, 8)  # single pdf
                                 
                                 usedmodel <- ''
                                 best_fits <- curve_fits[[endpoint]]
                                 all_experiments <- names(best_fits)
                                 for(experiment_id in all_experiments){
                                   
                                   visual_lines <- list()
                                   visual_points <- list()
                                   
                                   endpoint_bmc <- data.frame('bmr'=double(), 'lower_bmc'=double(), 'bmc'=double(), 'upper_bmc'=double())
                                   xAchse <- NULL
                                   
                                   best_fit <- best_fits[[experiment_id]]
                                   
                                   # prepare x-axis
                                   xlim_start = min(best_fit$origData[["dose"]])
                                   xlim_end <- max(best_fit$origData[["dose"]])*self$extrapolation
                                   
                                   xAchse <- exp(seq(log(xlim_start), log(xlim_end), length = self$grid_size))
                                   xAchse[1] <- xlim_start
                                   xAchse[self$grid_size] <- xlim_end
                                   
                                   best_fit_model <- best_fit$model
                                   
                                   if(is.na(best_fit$model)){ # protocol
                                     usedmodel <- 'no-fit'
                                     paste('no fit possible', experiment_id, self$compound, endpoint, sep = ' ') %push% self # log this process
                                   }
                                   else{
                                     
                                     usedmodel <- best_fit$model
                                     
                                     if(best_fit$model == 'flagged'){ # protocol
                                       paste('no plausible fit found for', experiment_id, self$compound, endpoint, sep = ' ') %push% self # log this process
                                     }
                                     
                                     neue_daten <- data.frame(
                                       dose = xAchse
                                     )
                                     
                                     plate <- predict(best_fit$drcObject, newdata =  neue_daten, interval = 'confidence')
                                     plate <- cbind(neue_daten, plate)
                                     colnames(plate) <- c('dose', 'best', 'lower', 'upper')
                                     
                                     
                                     dtpt_file_path <- paste0(self$output_folder, '/concentration-response-curves-data')
                                     dir.create(dtpt_file_path)
                                     
                                     dtpt_file_path <- paste0(dtpt_file_path, '/', self$compound)
                                     dir.create(dtpt_file_path)
                                     write.csv(plate, file = file.path(dtpt_file_path, paste0(experiment_id, '_', new_endpoint_pdf_name, '.csv')))
                                     
                                     
                                     # confidence bands and predicted fit provided as 
                                     # plate = data.frame(dose =, best =, lower =, upper = )
                                     
                                     supp_col <- data.frame('Endpoint' = new_endpoint_name, 'experiment' = experiment_id, 'flagged' = ifelse(best_fit_model == 'flagged', .2, 1))
                                     plate_data <- cbind(plate, supp_col)
                                     
                                     # plot
                                     visual_lines <- plate_data
                                     
                                     # bmc
                                     endpoint_bmc <- self$bmc_calc(for_endpoint=endpoint, plate=plate, bmr=bmr)
                                   }
                                   
                                   visual_points <- best_fit$origData
                                   visual_points$Endpoint <- new_endpoint_name
                                   
                                   # Dunnett Test
                                   dunnett.test.res <- dunnettSDP(visual_points)
                                   
                                   endpoint_bmc <- as.data.frame(endpoint_bmc)
                                   
                                   if(is.null(xAchse)){
                                     print('Hier geskippt.')
                                     next
                                   }
                                   
                                   xAchse[1] <- 0
                                   
                                   # add bmr data to visual line
                                   bmc_vertical_lower <- data.frame('dose'=double(), 'resp'=double(), 'Endpoint'=character(), 'flagged'=numeric())
                                   bmc_vertical <- data.frame('dose'=double(), 'resp'=double(), 'Endpoint'=character(), 'flagged'=numeric())
                                   bmc_vertical_upper <- data.frame('dose'=double(), 'resp'=double(), 'Endpoint'=character(), 'flagged'=numeric())
                                   
                                   for(indx in 1:length(bmr)){
                                     bmr_j <- bmr[[indx]]
                                     bmr_row <- data.frame('dose'= xAchse, 'best' = 100-bmr_j, 'lower' = NA, 'upper' = NA, 'Endpoint' = paste0('BMR', bmr_j), 'experiment' = NA, 'flagged' = 1)
                                     visual_lines <- rbind(visual_lines, bmr_row)
                                     
                                     if(nrow(endpoint_bmc) > 0){
                                       bmc_vertical_lower_row <- data.frame('dose'= endpoint_bmc[endpoint_bmc$"bmr" == 100-bmr_j, 'lower_bmc'], 'resp' = c(-Inf,100-bmr_j), 'Endpoint' = paste0('BMR', bmr_j),  'flagged' = 1)
                                       bmc_vertical_lower <- rbind(bmc_vertical_lower, bmc_vertical_lower_row)
                                       bmc_vertical_row <- data.frame('dose'= endpoint_bmc[endpoint_bmc$"bmr" == 100-bmr_j, 'bmc'][1], 'resp' = c(-Inf,100-bmr_j), 'Endpoint' = paste0('BMR', bmr_j),  'flagged' = 1)
                                       bmc_vertical <- rbind(bmc_vertical, bmc_vertical_row)
                                       bmc_vertical_upper_row <- data.frame('dose'= endpoint_bmc[endpoint_bmc$"bmr" == 100-bmr_j, 'upper_bmc'], 'resp' = c(-Inf,100-bmr_j), 'Endpoint' = paste0('BMR', bmr_j),  'flagged' = 1)
                                       bmc_vertical_upper <- rbind(bmc_vertical_upper, bmc_vertical_upper_row)
                                     }
                                   }
                                   
                                   colnames(bmc_vertical_lower) <- c('dose', 'resp', 'Endpoint', 'flagged')
                                   colnames(bmc_vertical) <- c('dose', 'resp', 'Endpoint', 'flagged')
                                   colnames(bmc_vertical_upper) <- c('dose', 'resp', 'Endpoint', 'flagged')
                                   
                                   # plot curves
                                   p <- ggplot2::ggplot() + ggplot2::ylim(0, self$yaxis_lim) + ggplot2::scale_x_log10() +
                                     ggplot2::ylab(as.expression(bquote(bold("% of control")))) +
                                     ggplot2::xlab(as.expression(bquote(bold(.(self$compound) ~~ "["~mu~M~"]") ) )) + 
                                     ggplot2::labs(title = experiment_id, caption = usedmodel) +
                                     ggplot2::theme_classic(base_size=16) +
                                     ggplot2::theme(legend.position = 'top',
                                                    legend.direction='vertical',
                                                    legend.title = ggplot2::element_blank(),
                                                    legend.text = ggplot2::element_text(face = 'bold'),
                                                    axis.line = ggplot2::element_line(size = 1, linetype = "solid"),
                                                    axis.text = ggplot2::element_text(face = 'bold'))
                                   
                                   logicvec <- new_endpoints_names %in% new_endpoint_name
                                   logicvec <- c(logicvec, rep(TRUE, length(bmr)))
                                   
                                   if(!all(is.na(visual_lines['lower']))){
                                     p <- p + ggplot2::geom_ribbon(data = visual_lines, mapping = ggplot2::aes(x=dose, y=best, group=Endpoint, ymin=lower, ymax=upper), fill = 'yellow')
                                   }
                                   
                                   p <- p + ggplot2::geom_line(data = visual_lines, mapping = ggplot2::aes(x=dose, y=best, group=Endpoint, color=Endpoint, linetype=Endpoint, alpha=Endpoint, size=Endpoint)) +
                                     ggplot2::scale_color_manual(values=group_colors[logicvec]) +
                                     ggplot2::scale_linetype_manual(values=group_linetype[logicvec]) + 
                                     ggplot2::scale_alpha_manual(values=group_alpha[logicvec]) +
                                     ggplot2::scale_size_manual(values=group_size[logicvec])
                                   
                                   # origData
                                   if(experiments == 'merged'){
                                     p <- p + ggplot2::geom_point(data = visual_points, mapping = ggplot2::aes(x=dose, y=endpoint, shape=experiment), size=4, color='#800080')
                                     
                                     if(nrow(dunnett.test.res) > 0){
                                       dunnett.test.res[, ('y'):= 150.0]
                                       p <- p + ggplot2::geom_text(data = dunnett.test.res, mapping = ggplot2::aes(x=dose, y=y, label = p.signif), size=10, face = 'bold')
                                     }
                                     
                                     
                                   }
                                   else{
                                     p <- p + ggplot2::geom_point(data = visual_points, mapping = ggplot2::aes(x=dose, y=endpoint), shape=20, size=4, color='#800080')
                                   }
                                   
                                   
                                   dtpt_file_path <- paste0(self$output_folder, '/concentration-response-points')  
                                   dir.create(dtpt_file_path)
                                   
                                   dtpt_file_path <- paste0(dtpt_file_path, '/', self$compound)
                                   dir.create(dtpt_file_path)
                                   
                                   tmp_visual_points <- convert_table(method = experiments, visual_points)
                                   write.csv(tmp_visual_points, file = file.path(dtpt_file_path, paste0(experiment_id, '_', new_endpoint_pdf_name, '.csv')))
                                   
                                   # bmc
                                   if(nrow(endpoint_bmc) > 0){
                                     p <- p + ggplot2::geom_line(data = bmc_vertical_lower, mapping = ggplot2::aes(x=dose, y=resp, group=Endpoint, color=Endpoint, linetype=Endpoint, alpha=Endpoint, size=Endpoint)) 
                                     p <- p + ggplot2::geom_line(data = bmc_vertical, mapping = ggplot2::aes(x=dose, y=resp, group=Endpoint, color=Endpoint, linetype=Endpoint, alpha=Endpoint, size=Endpoint)) 
                                     p <- p + ggplot2::geom_line(data = bmc_vertical_upper, mapping = ggplot2::aes(x=dose, y=resp, group=Endpoint, color=Endpoint, linetype=Endpoint, alpha=Endpoint, size=Endpoint)) 
                                   }
                                   
                                   visual_points <- Rmisc::summarySE(visual_points, measurevar="endpoint", groupvars=c("Endpoint","dose"))
                                   
                                   p <- p + ggplot2::geom_errorbar(data = visual_points, mapping = ggplot2::aes(x=dose, y=endpoint, ymin=endpoint-se, ymax=endpoint+se), width=.05, size=.5, alpha=.7)
                                   
                                   p <- p + ggplot2::geom_point(data = visual_points, mapping = ggplot2::aes(x=dose, y=endpoint), alpha=.7)
                                   
                                   p <- p + ggplot2::annotation_logticks(sides = "b", short = ggplot2::unit(0.1, "cm"), mid = ggplot2::unit(0.1, "cm"), long = ggplot2::unit(0.1, "cm"))
                                   plot(p)
                                 }
                                 
                                 
                                 dev.off()
                               }
                               
                               
                             },
                             
                             
                             #' @description
                             #' Print curves of the endpoints.
                             #' @param experiments merged or single option.
                             #' @param for_endpoint character or list of endpoints. If not restricted, for all endpoints will be printed.
                             #' @param normal TRUE (default) or FALSE for renormalisation of the data regarding the smallest dose.
                             #' @param bmr double or list of double with the benchmark-response
                             #' @return Plot of the data regarding to the endpoints (and experiments) with the benchmark-response and benchmark-dose.                         
                             classification_plot = function(experiments = c('merged'), for_endpoint = NULL, specific_vs_unspecific_table = NA, classification_rules_table = NA){
                               # parameter
                               experiments <- match.arg(experiments)
                               
                               if((typeof(for_endpoint) == 'list') && (nrow(for_endpoint) > 0)){ # no slashes allowed in pdfname column
                                 endpoints_final_names_table <- for_endpoint
                               }
                               else{
                                 endpoints_final_names_table <- data.table('Name of PDF file' = gsub("/", "", self$get_endpoints()), 'Name in PDF'=self$get_endpoints(), 'Input'=self$get_endpoints())
                               }
                               
                               if((typeof(specific_vs_unspecific_table) != 'list') || nrow(specific_vs_unspecific_table) == 0) return('Specific vs Unspecific declaring reference table is missing! No classification could be applied.');
                               
                               if((typeof(classification_rules_table) != 'list') || nrow(classification_rules_table) == 0) return('Classification Rules Table table is missing! No classification could be applied.');
                               
                               N <- nrow(endpoints_final_names_table)
                               
                               result <- data.table() # for classification
                               classification_row <- list()
                               
                               path <- paste0(self$output_folder, '/classification-crcs')
                               dir.create(path)
                               path <- paste0(path, '/', self$compound)
                               dir.create(path)
                               
                               # assumption correct/consistent mode (variable 'experiment')
                               curve_fits <- self$crc_dt
                               
                               endpoints <- self$get_endpoints()
                               specific_endpoints <- specific_vs_unspecific_table[, get('Specific endpoint')]
                               
                               for(i in 1:N){
                                 endpoint <- endpoints_final_names_table[i, get('Input')]
                                 
                                 if(!(endpoint %in% endpoints)) next;
                                 if(!(endpoint %in% specific_endpoints)) next;
                                 
                                 new_endpoint_name <- endpoints_final_names_table[i, get('Name in PDF')]
                                 new_endpoint_pdf_name <- endpoints_final_names_table[i, get('Name of PDF file')]
                                 new_endpoint_pdf_name <- gsub("/", " ", new_endpoint_pdf_name) # to ensure valid file name of the output pdf
                                 
                                 best_fits <- curve_fits[[endpoint]]
                                 best_fit <- best_fits[['merged']]
                                 usedmodel <- ''
                                 best_fit_model <- best_fit$model
                                 
                                 # +++  +++ +++  +++ +++  +++
                                 # prepare xaxis of common plot. specified by xaxis interval of (specific) endpoint
                                 xlim_start = min(best_fit$origData[["dose"]])
                                 xlim_end <- max(best_fit$origData[["dose"]])*self$extrapolation
                                 xAchse <- exp(seq(log(xlim_start), log(xlim_end), length = self$grid_size))
                                 xAchse[1] <- xlim_start
                                 xAchse[self$grid_size] <- xlim_end
                                 neue_daten <- data.frame(
                                   dose = xAchse
                                 )
                                 # +++  +++ +++  +++ +++  +++
                                 
                                 
                                 
                                 classification_row[['Specific endpoint']] <- endpoint
                                 
                                 visual_lines_specific <- data.frame()
                                 visual_points_specific <- data.frame()
                                 
                                 # +++  +++ +++  +++ +++  +++
                                 # specific endpoint plot data
                                 if(is.na(best_fit$model)){
                                   
                                   usedmodel <- 'no-fit'
                                   classification_row[['Trend specific']] <- 'NA'
                                   
                                 }
                                 else{
                                   
                                   usedmodel <- best_fit$model
                                   plate <- predict(best_fit$drcObject, newdata =  neue_daten, interval = 'confidence')
                                   plate <- cbind(neue_daten, plate)
                                   colnames(plate) <- c('dose', 'best', 'lower', 'upper')
                                   
                                   plate_specific <- plate
                                   
                                   supp_col <- data.frame('Endpoint' = new_endpoint_name, 'experiment' = 'merged', 'flagged' = ifelse(best_fit_model == 'flagged', .2, 1))
                                   plate_data <- cbind(plate, supp_col)
                                   
                                   # plot
                                   visual_lines_specific <- plate_data
                                   
                                   
                                   if( (plate[1, 'best'] - plate[self$grid_size, 'best']) > 0 ){
                                     classification_row[['Trend specific']] <- 'reduction'
                                   }
                                   else if( (plate[1, 'best'] - plate[self$grid_size, 'best']) < 0 ){
                                     classification_row[['Trend specific']] <- 'induction'
                                   }
                                   else{
                                     classification_row[['Trend specific']] <- 'no-effect'
                                   }
                                   
                                 }
                                 
                                 usedmodel_specific <- usedmodel
                                 visual_points_specific <- best_fit$origData
                                 visual_points_specific$Endpoint <- new_endpoint_name
                                 
                                 dunnett.test.res <- dunnettSDP(visual_points_specific)
                                 if(nrow(dunnett.test.res) > 0) dunnett.test.res <- dunnett.test.res[!is.na(get('q-values'))]; # na.omit
                                 
                                 # print(paste0('specific ', endpoint))
                                 # View(dunnett.test.res)
                                 
                                 # Specific endpoint viability related?
                                 classification_row[['Specific endpoint viability related']] <- specific_vs_unspecific_table[get('Specific endpoint') == endpoint, get('Specific endpoint viability related')][1]
                                 
                                 
                                 # Specific endpoint highest concentration is significant?
                                 if(nrow(dunnett.test.res) > 0){
                                   dunndose <- dunnett.test.res[, get('dose')]
                                   dunndose <- max(dunndose)
                                   q_val_max_dose <- dunnett.test.res[get('dose') == dunndose, get('q-values')][1]
                                   ifelse(q_val_max_dose <= 0.05, classification_row[['Specific endpoint highest concentration is significant']] <- 'yes', classification_row[['Specific endpoint highest concentration is significant']] <- 'no')
                                 }
                                 else{
                                   classification_row[['Specific endpoint highest concentration is significant']] <- 'NA'
                                 }
                                 
                                 
                                 
                                 # +++  +++ +++  +++ +++  +++
                                 # visual_points and visual_lines must be extended by unspecific endpoints
                                 
                                 
                                 # +++  +++ +++  +++ +++  +++
                                 table_focus <- specific_vs_unspecific_table[get('Specific endpoint') == endpoint]
                                 M <- nrow(table_focus)
                                 
                                 
                                 for(j in 1:M){
                                   unspecific_endpoint <- table_focus[j, get('Unspecific endpoint')]
                                   
                                   if(!(unspecific_endpoint %in% endpoints)) next;
                                   if(!(unspecific_endpoint %in% endpoints_final_names_table[, get('Input')])) next;
                                   
                                   classification_row[['Unspecific endpoint']] <- unspecific_endpoint
                                   
                                   endpoint_bmc <- data.frame('bmr'=double(), 'lower_bmc'=double(), 'bmc'=double(), 'upper_bmc'=double())
                                   unspecific_endpoint_bmc <- data.frame('bmr'=double(), 'lower_bmc'=double(), 'bmc'=double(), 'upper_bmc'=double())
                                   
                                   # specific bmc
                                   specific_bmr <- strsplit(table_focus[j, get('Specific BMR')], ', ')[[1]]
                                   specific_bmr <- as.list(as.double(specific_bmr))
                                   if(usedmodel_specific != 'no-fit'){
                                     endpoint_bmc <- self$bmc_calc(for_endpoint=endpoint, plate=plate_specific, bmr=specific_bmr)
                                   }
                                   else{
                                     endpoint_bmc <- data.frame('endpoint'=endpoint, 'bmr'=as.double(specific_bmr), 'lower_bmc'=as.numeric(NA), 'bmc'=as.numeric(NA), 'upper_bmc'=as.numeric(NA))
                                   }
                                   
                                   unspecific_bmr <- strsplit(table_focus[j, get('Unspecific BMR')], ', ')[[1]]
                                   unspecific_bmr <- as.list(as.double(unspecific_bmr))
                                   
                                   best_fits <- curve_fits[[unspecific_endpoint]]
                                   best_fit <- best_fits[['merged']]
                                   usedmodel <- ''
                                   best_fit_model <- best_fit$model
                                   
                                   visual_lines_unspecific <- data.frame()
                                   visual_points_unspecific <- data.frame()
                                   
                                   if(is.na(best_fit$model)){
                                     
                                     usedmodel <- 'no-fit'
                                     classification_row[['Trend unspecific']] <- 'NA'
                                     
                                     unspecific_endpoint_bmc <- data.frame('endpoint'=unspecific_endpoint, 'bmr'=as.double(unspecific_bmr), 'lower_bmc'=as.numeric(NA), 'bmc'=as.numeric(NA), 'upper_bmc'=as.numeric(NA))
                                   }
                                   else{
                                     
                                     usedmodel <- best_fit$model
                                     plate <- predict(best_fit$drcObject, newdata =  neue_daten, interval = 'confidence')
                                     plate <- cbind(neue_daten, plate)
                                     colnames(plate) <- c('dose', 'best', 'lower', 'upper')
                                     
                                     supp_col <- data.frame('Endpoint' = unspecific_endpoint, 'experiment' = 'merged', 'flagged' = ifelse(best_fit_model == 'flagged', .2, 1))
                                     plate_data <- cbind(plate, supp_col)
                                     
                                     # plot
                                     visual_lines_unspecific <- plate_data
                                     
                                     # trend of fit
                                     if( (plate[1, 'best'] - plate[self$grid_size, 'best']) > 0 ){
                                       classification_row[['Trend unspecific']] <- 'reduction'
                                     }
                                     else if( (plate[1, 'best'] - plate[self$grid_size, 'best']) < 0 ){
                                       classification_row[['Trend unspecific']] <- 'induction'
                                     }
                                     else{
                                       classification_row[['Trend unspecific']] <- 'no-effect'
                                     }
                                     
                                     # bmc
                                     unspecific_endpoint_bmc <- self$bmc_calc(for_endpoint=unspecific_endpoint, plate=plate, bmr=unspecific_bmr)
                                   }
                                   
                                   visual_points_unspecific <- best_fit$origData
                                   visual_points_unspecific$Endpoint <- unspecific_endpoint
                                   
                                   endpoint_bmc <- as.data.frame(endpoint_bmc)
                                   unspecific_endpoint_bmc <- as.data.frame(unspecific_endpoint_bmc)
                                   
                                   # merge specific & unspecific
                                   visual_lines <- rbind(visual_lines_specific, visual_lines_unspecific)
                                   visual_points <- rbind(visual_points_specific, visual_points_unspecific)
                                   
                                   bmr_types <- list(list('bmrs' = specific_bmr, 'bmcs' = endpoint_bmc, 'type' = 'specific'),
                                                     list('bmrs' = unspecific_bmr, 'bmcs' = unspecific_endpoint_bmc, 'type' = 'unspecific')
                                   )
                                   
                                   
                                   N_bmrs <- length(bmr_types[[1]][['bmrs']])
                                   
                                   for(k in 1:N_bmrs){
                                     # PDF for each row
                                     grDevices::pdf(file = file.path(path, paste0(new_endpoint_pdf_name,'_vs_', unspecific_endpoint, k, '.pdf')), 9, 8)  # single pdf
                                     
                                     # classification 
                                     classification_row_further <- list()
                                     classification_row_further[['Specific endpoint']] <- classification_row[['Specific endpoint']]
                                     classification_row_further[['Unspecific endpoint']] <- classification_row[['Unspecific endpoint']]
                                     classification_row_further[['Trend specific']] <- classification_row[['Trend specific']]
                                     classification_row_further[['Specific endpoint viability related']] <- classification_row[['Specific endpoint viability related']]
                                     classification_row_further[['Specific endpoint highest concentration is significant']] <- classification_row[['Specific endpoint highest concentration is significant']]
                                     classification_row_further[['Trend unspecific']] <- classification_row[['Trend unspecific']]
                                     
                                     
                                     bmr_groups <- c()
                                     bmc_groups <- c()
                                     bmc_verwalter <- list()
                                     # add bmr data to visual line
                                     bmc_vertical_lower <- data.frame('dose'=double(), 'resp'=double(), 'Endpoint'=character(), 'flagged'=numeric())
                                     bmc_vertical <- data.frame('dose'=double(), 'resp'=double(), 'Endpoint'=character(), 'flagged'=numeric())
                                     bmc_vertical_upper <- data.frame('dose'=double(), 'resp'=double(), 'Endpoint'=character(), 'flagged'=numeric())
                                     
                                     for(l in 1:length(bmr_types)){
                                       ctype <- bmr_types[[l]][['type']]
                                       bmr <- bmr_types[[l]][['bmrs']]
                                       cbmcs <- bmr_types[[l]][['bmcs']]
                                       
                                       bmr_j <- bmr[[k]]
                                       bmr_row <- data.frame('dose'= xAchse, 'best' = 100-bmr_j, 'lower' = NA, 'upper' = NA, 'Endpoint' = paste0('BMR', bmr_j, ctype), 'experiment' = NA, 'flagged' = 1)
                                       visual_lines <- rbind(visual_lines, bmr_row)
                                       if(!(paste0('BMR', bmr_j, ctype) %in% bmr_groups)) bmr_groups <- c(bmr_groups, paste0('BMR', bmr_j, ctype));
                                       
                                       bmc_verwalter[[ctype]] <- list('lbmc'=NA, 'bmc'=NA, 'ubmc'=NA)
                                       
                                       if(nrow(cbmcs) > 0){
                                         # lower
                                         lower.bmc <- cbmcs[cbmcs$"bmr" == 100-bmr_j, 'lower_bmc']
                                         for(current_bmc in lower.bmc){
                                           bmc_vertical_lower_row <- data.frame('dose'= current_bmc, 'resp' = c(-Inf,100-bmr_j), 'Endpoint' = paste0('BMC', bmr_j, ctype),  'flagged' = 1)
                                           bmc_vertical_lower <- rbind(bmc_vertical_lower, bmc_vertical_lower_row)
                                           
                                           lower.bmc.interval <- current_bmc
                                         }
                                         
                                         # 
                                         estimate.bmc <- cbmcs[cbmcs$"bmr" == 100-bmr_j, 'bmc']
                                         estimate.bmc.interval <- NA
                                         for(current_bmc in estimate.bmc){
                                           bmc_vertical_row <- data.frame('dose'= current_bmc, 'resp' = c(-Inf,100-bmr_j), 'Endpoint' = paste0('BMC', bmr_j, ctype),  'flagged' = 1)
                                           bmc_vertical <- rbind(bmc_vertical, bmc_vertical_row)
                                           
                                           estimate.bmc.interval <- current_bmc
                                         }
                                         
                                         # upper
                                         upper.bmc <- cbmcs[cbmcs$"bmr" == 100-bmr_j, 'upper_bmc']
                                         for(current_bmc in upper.bmc){
                                           bmc_vertical_upper_row <- data.frame('dose'= current_bmc, 'resp' = c(-Inf,100-bmr_j), 'Endpoint' = paste0('BMC', bmr_j, ctype),  'flagged' = 1)
                                           bmc_vertical_upper <- rbind(bmc_vertical_upper, bmc_vertical_upper_row)
                                           
                                           upper.bmc.interval <- current_bmc
                                         }
                                         
                                         if(!(paste0('BMC', bmr_j, ctype) %in% bmr_groups)) bmc_groups <- c(bmc_groups, paste0('BMC', bmr_j, ctype));
                                         
                                         # bmc-interval
                                         if(!is.na(estimate.bmc.interval)){
                                           classification_row_further[[paste('BMC', ctype, 'available', sep = ' ')]] <- 'yes'
                                           
                                           # ifelse(is.na(lower.bmc.interval), lbmc <- estimate.bmc.interval, lbmc <- lower.bmc.interval)
                                           ifelse(is.na(lower.bmc.interval), lbmc <- xAchse[1], lbmc <- lower.bmc.interval)
                                           ifelse(is.na(upper.bmc.interval), ubmc <- xAchse[self$grid_size], ubmc <- upper.bmc.interval)
                                           bmc_interval <- data.frame('dose'= xAchse[(lbmc <= xAchse) & (ubmc >= xAchse)], 'best' = 0, 'lower' = NA, 'upper' = NA, 'Endpoint' = paste0('BMC', bmr_j, ctype), 'experiment' = NA, 'flagged' = 1)
                                           
                                           bmc_verwalter[[ctype]] <- list('lbmc'=lbmc, 'bmc'=estimate.bmc.interval, 'ubmc'=ubmc)
                                           
                                           visual_lines <- rbind(visual_lines, bmc_interval)
                                         }
                                         else{
                                           classification_row_further[[paste('BMC', ctype, 'available', sep = ' ')]] <- 'no'
                                         }
                                         
                                       }
                                       
                                     }
                                     
                                     colnames(bmc_vertical_lower) <- c('dose', 'resp', 'Endpoint', 'flagged')
                                     colnames(bmc_vertical) <- c('dose', 'resp', 'Endpoint', 'flagged')
                                     colnames(bmc_vertical_upper) <- c('dose', 'resp', 'Endpoint', 'flagged')
                                     
                                     # +++  +++ +++  +++ +++  +++
                                     classification_row_further[['BMR specific']] <- bmr_types[[1]][['bmrs']][[k]]
                                     classification_row_further[['BMR unspecific']] <- bmr_types[[2]][['bmrs']][[k]]
                                     
                                     classification_row_further[['BMC specific']] <- bmc_verwalter[['specific']][['bmc']]
                                     classification_row_further[['BMCL specific']] <- bmc_verwalter[['specific']][['lbmc']]
                                     classification_row_further[['BMCU specific']] <- bmc_verwalter[['specific']][['ubmc']]
                                     classification_row_further[['BMC unspecific']] <- bmc_verwalter[['unspecific']][['bmc']]
                                     classification_row_further[['BMCL unspecific']] <- bmc_verwalter[['unspecific']][['lbmc']]
                                     classification_row_further[['BMCU unspecific']] <- bmc_verwalter[['unspecific']][['ubmc']]
                                     
                                     
                                     # BMC specific upper limit within test range ?
                                     if(!is.na(bmc_verwalter[['specific']][['bmc']])){
                                       if(bmc_verwalter[['specific']][['ubmc']] <= xlim_end/self$extrapolation){
                                         classification_row_further[['BMC specific upper limit within test range']] <- 'yes'
                                       }
                                       else{
                                         classification_row_further[['BMC specific upper limit within test range']] <- 'no'
                                       }
                                     }
                                     
                                     
                                     
                                     # Some further data for classification
                                     if( (!is.na(bmc_verwalter[['specific']][['bmc']])) &&  (!is.na(bmc_verwalter[['unspecific']][['bmc']])) ){
                                       
                                       # BMC unspecific > BMC specific?
                                       if( bmc_verwalter[['unspecific']][['bmc']] > bmc_verwalter[['specific']][['bmc']] ){
                                         classification_row_further[['BMC unspecific > BMC specific']] <- 'yes'
                                       }
                                       else{
                                         classification_row_further[['BMC unspecific > BMC specific']] <- 'no'
                                       }
                                       
                                       
                                       # Overlap-Ratio
                                       classification_row_further[['Overlap-Ratio']] <- (bmc_verwalter[['unspecific']][['lbmc']] - bmc_verwalter[['specific']][['ubmc']]) / (bmc_verwalter[['specific']][['ubmc']] - bmc_verwalter[['specific']][['lbmc']])
                                       
                                     }
                                     else{
                                       
                                       classification_row_further[['BMC unspecific > BMC specific']] <- 'NA'
                                       
                                       # Overlap-Ratio
                                       classification_row_further[['Overlap-Ratio']] <- 'NA'
                                       
                                     }
                                     
                                     
                                     
                                     
                                     
                                     # +++  +++ +++  +++ +++  +++
                                     
                                     
                                     
                                     # +++  +++ +++  +++ +++  +++
                                     # PLOT 
                                     p <- ggplot2::ggplot() + ggplot2::ylim(0, self$yaxis_lim) + ggplot2::scale_x_log10() +
                                       ggplot2::ylab(as.expression(bquote(bold("% of control")))) +
                                       ggplot2::xlab(as.expression(bquote(bold(.(self$compound) ~~ "["~mu~M~"]") ) )) + 
                                       ggplot2::labs(title = 'classification', caption = 'specific vs unspecific') +
                                       ggplot2::theme_classic(base_size=16) +
                                       ggplot2::theme(legend.position = 'top',
                                                      legend.direction='vertical',
                                                      legend.title = ggplot2::element_blank(),
                                                      legend.text = ggplot2::element_text(face = 'bold'),
                                                      axis.line = ggplot2::element_line(size = 1, linetype = "solid"),
                                                      axis.text = ggplot2::element_text(face = 'bold'))
                                     
                                     group_colors <- c('#00167B', '#bd9564', rep('#B9B7BD', length(bmr_groups)), '#7c8cd6', '#FFFF00')
                                     group_shapes <- c(18, 18, rep(NA, length(bmr_groups)), NA, NA)
                                     group_linetype <- c('solid', 'solid', 'dashed', 'dotdash', 'solid', 'solid')
                                     group_alpha <- c(1, 1, rep(.2, length(bmr_groups)), .6, .6)
                                     group_size <- c(1.6, 1.6, rep(.5, length(bmr_groups)), 1.6, 1.6)
                                     group_fill <- c('#7c8cd6', '#FFFF00', rep('#f8f8ff', length(bmr_groups)), rep('#f8f8ff', length(bmc_groups)))
                                     
                                     names(group_colors) <- c(new_endpoint_name, unspecific_endpoint, bmr_groups, bmc_groups)
                                     names(group_shapes) <- c(new_endpoint_name, unspecific_endpoint, bmr_groups, bmc_groups)
                                     names(group_linetype) <- c(new_endpoint_name, unspecific_endpoint, bmr_groups, bmc_groups)
                                     names(group_alpha) <- c(new_endpoint_name, unspecific_endpoint, bmr_groups, bmc_groups)
                                     names(group_size) <- c(new_endpoint_name, unspecific_endpoint, bmr_groups, bmc_groups)
                                     names(group_fill) <- c(new_endpoint_name, unspecific_endpoint, bmr_groups, bmc_groups)
                                     
                                     if(!all(is.na(visual_lines['lower']))){
                                       p <- p + ggplot2::geom_ribbon(data = visual_lines, mapping = ggplot2::aes(x=dose, y=best, group=Endpoint, ymin=lower, ymax=upper, fill = Endpoint), alpha=.5) +
                                         ggplot2::scale_fill_manual(values=group_fill)
                                     }
                                     
                                     p <- p + ggplot2::geom_line(data = visual_lines, mapping = ggplot2::aes(x=dose, y=best, group=Endpoint, color=Endpoint, linetype=Endpoint, alpha=Endpoint, size=Endpoint)) +
                                       ggplot2::scale_color_manual(values=group_colors) +
                                       ggplot2::scale_linetype_manual(values=group_linetype) + 
                                       ggplot2::scale_alpha_manual(values=group_alpha) +
                                       ggplot2::scale_size_manual(values=group_size)
                                     
                                     # origData
                                     if(experiments == 'merged'){
                                       p <- p + ggplot2::geom_point(data = visual_points, mapping = ggplot2::aes(x=dose, y=endpoint, group=Endpoint, shape=experiment, color=Endpoint), size=4)
                                       
                                       if(nrow(dunnett.test.res) > 0){
                                         dunnett.test.res[, ('y'):= 150.0]
                                         p <- p + ggplot2::geom_text(data = dunnett.test.res, mapping = ggplot2::aes(x=dose, y=y, label = p.signif), size=10, face = 'bold')
                                       }
                                       
                                     }
                                     
                                     # bmc
                                     if(nrow(bmc_vertical) > 0){
                                       p <- p + ggplot2::geom_line(data = bmc_vertical_lower, mapping = ggplot2::aes(x=dose, y=resp, group=Endpoint, color=Endpoint, linetype=Endpoint, alpha=Endpoint, size=Endpoint)) 
                                       # p <- p + geom_line(data = bmc_vertical, mapping = aes(x=dose, y=resp, group=Endpoint, color=Endpoint, linetype=Endpoint, alpha=Endpoint, size=Endpoint)) 
                                       p <- p + ggplot2::geom_line(data = bmc_vertical_upper, mapping = ggplot2::aes(x=dose, y=resp, group=Endpoint, color=Endpoint, linetype=Endpoint, alpha=Endpoint, size=Endpoint)) 
                                     }
                                     
                                     
                                     if(FALSE){
                                       visual_points_se <- Rmisc::summarySE(visual_points, measurevar="endpoint", groupvars=c("Endpoint","dose"))
                                       p <- p + ggplot2::geom_errorbar(data = visual_points_se, mapping = ggplot2::aes(x=dose, y=endpoint, ymin=endpoint-se, ymax=endpoint+se, group=Endpoint, color=Endpoint), width=.05, size=.5, alpha=.7)
                                       p <- p + ggplot2::geom_point(data = visual_points_se, mapping = ggplot2::aes(x=dose, y=endpoint, group=Endpoint, color=Endpoint), alpha=.7)
                                       
                                     }
                                     
                                     p <- p + ggplot2::annotation_logticks(sides = "b", short = ggplot2::unit(0.1, "cm"), mid = ggplot2::unit(0.1, "cm"), long = ggplot2::unit(0.1, "cm"))
                                     
                                     plot(p)
                                     dev.off()
                                     
                                     # +++  +++ +++  +++ +++  +++
                                     # print(colnames(as.data.frame(classification_row_further)))
                                     setDT(classification_row_further)
                                     result <- rbind(result, classification_row_further, fill = TRUE)
                                     # +++  +++ +++  +++ +++  +++
                                     
                                     
                                   }
                                   
                                 }
                                 
                               }
                               
                               
                               # apply classification on result according 'classification'
                               result[, ('Classification'):= '']
                               result[, ('Comment'):= '']
                               
                               # View(result)
                               N <- nrow(result)
                               M <- nrow(classification_rules_table)
                               rule.attrs <- colnames(classification_rules_table)
                               for(i in 1:N){
                                 for(j in 1:M){
                                   crule <- classification_rules_table[j, ]
                                   
                                   # determine rules of interest
                                   attrs_of_interest <- c()
                                   for(rule.attr in rule.attrs){
                                     if(rule.attr %in% c('Overlap-Ratio Interval', 'Classification', 'Comment')) next;
                                     
                                     if(is.na(crule[, get(rule.attr)]) || (crule[, get(rule.attr)] == '')) next;
                                     attrs_of_interest <- c(attrs_of_interest, rule.attr)
                                   }
                                   
                                   # print(attrs_of_interest)
                                   if(!all(attrs_of_interest %in% colnames(result[i, ]))) next;
                                   
                                   # View(result[i, ..attrs_of_interest])
                                   # View(crule[, ..attrs_of_interest])
                                   
                                   if( !all( as.character(result[i, ..attrs_of_interest]) == as.character(crule[, ..attrs_of_interest]) ) ) next;
                                   
                                   # now check Overlap-Ratio Interval
                                   if((!is.na(crule[, get('Overlap-Ratio Interval')])) && (crule[, get('Overlap-Ratio Interval')] != '')){
                                     
                                     if(is.nan(result[i, get('Overlap-Ratio')]) || is.na(result[i, get('Overlap-Ratio')])) next;
                                     
                                     cinterval <- crule[, get('Overlap-Ratio Interval')]
                                     first_bracket <- substr(cinterval, 1, 1)
                                     last_bracket <- substr(cinterval, nchar(cinterval), nchar(cinterval))
                                     
                                     cinterval_wo_brackets <- gsub(pattern = "\\(|\\[|\\)|\\]", replacement = "", cinterval)
                                     cinterval_vec <- as.numeric(strsplit(cinterval_wo_brackets, ', ')[[1]])
                                     
                                     if(first_bracket == '('){
                                       
                                       if(as.numeric(result[i, get('Overlap-Ratio')]) <= cinterval_vec[1]) next;
                                       
                                     }
                                     else if(first_bracket == '['){
                                       
                                       if(as.numeric(result[i, get('Overlap-Ratio')]) < cinterval_vec[1]) next;
                                       
                                     }
                                     else{  # syntax wrong in table
                                       
                                       next
                                       
                                     }
                                     
                                     
                                     
                                     if(last_bracket == ')'){
                                       
                                       if(as.numeric(result[i, get('Overlap-Ratio')]) >= cinterval_vec[2]) next;
                                       
                                     }
                                     else if(last_bracket == ']'){
                                       
                                       if(as.numeric(result[i, get('Overlap-Ratio')]) > cinterval_vec[2]) next;
                                       
                                     }
                                     else{  # syntax wrong in table
                                       
                                       next
                                       
                                     }
                                     
                                   }
                                   
                                   classification_result <- crule[, get('Classification')]
                                   comment_result <- crule[, get('Comment')]
                                   result[i, ('Classification'):= classification_result]
                                   result[i, ('Comment'):= comment_result]
                                   break
                                   
                                 }
                                 
                                 
                               }
                               
                               # save
                               path <- paste0(self$output_folder, '/classification')
                               dir.create(path)
                               write.table(result, file = paste0(path, '/', self$compound, '.csv'), sep=';', row.names = FALSE)
                               
                             },
                             
                             #' @description
                             #' Starts a full run.
                             #' @param task xlsx file containing the tasks and parameter.
                             #' @return 1                         
                             run = function(task_file = NA){
                               
                               # load task_file
                               tasks <- tryCatch({
                                 tmp <- readxl::read_excel(file.path(task_file), sheet = "Task")
                                 tmp
                               }, 
                               error=function(cond) {
                                 print('Error reading task file.' )
                                 return(data.table())
                               },
                               warning=function(cond) {
                                 print(cond$message) 
                                 return(data.table())
                               }
                               )
                               
                               tasks <- as.data.frame(tasks)
                               tasks <- tasks[!is.na(tasks$Task),]
                               setDT(tasks)
                               utils::write.table(tasks, file = paste0(self$output_folder, '/Tasks.csv'), sep=';', row.names = FALSE)
                               tasks <- tasks[get('Apply') == 'yes']
                               
                               N_tasks <- nrow(tasks)
                               for(i in 1:N_tasks){
                                 ctask <- tasks[i, get('Task')] 
                                 switch (trim(ctask),
                                         "compute quotient parameter" = {
                                           quotient_table <- readxl::read_excel(file.path(task_file), sheet = "Quotients")
                                           setDT(quotient_table)
                                           
                                           self$add_quotient_parameter(quotients = quotient_table)
                                         },
                                         "custom cutoff" = {
                                           cutoff_table <- readxl::read_excel(file.path(task_file), sheet = "CustomCutOff")
                                           setDT(cutoff_table)
                                           
                                           self$custom_cut_off(cutoff_table = cutoff_table)
                                         },
                                         "ensure minimal data requirements" = {
                                           self$control_replicates()
                                           self$conditions()
                                         },
                                         "background normalization" = {
                                           bgt <- readxl::read_excel(file.path(task_file), sheet = "BackgroundNormalization")
                                           setDT(bgt)
                                           
                                           self$background_subtraction(background_table = bgt)
                                         },
                                         "control- + dynamic range normalization" = {
                                           dr_table <- readxl::read_excel(file.path(task_file), sheet = "DynamicRangeNormalization")
                                           setDT(dr_table)
                                           
                                           self$control_normalization(dr_endpoints_table = dr_table) 
                                         },
                                         "only control normalization" = {
                                           self$control_normalization() 
                                         },
                                         "fit data + concentration response curves + compute bmc estimates" = {
                                           bmr_file <- readxl::read_excel(file.path(task_file), sheet = "BMRs")
                                           setDT(bmr_file)
                                           
                                           endpoints_table <- readxl::read_excel(file.path(task_file), sheet = "Endpoints")
                                           setDT(endpoints_table)
                                           
                                           # self$fit_dt(experiments = 'merged', endpoints_table = endpoints_table, rnrm = 'yes')
                                           self$fit_dt(experiments = 'merged', endpoints_table = endpoints_table)
                                           self$plotc(experiments = 'merged', for_endpoint = endpoints_table, bmr = bmr_file)
                                         },
                                         "fit data + concentration response curves + compute bmc estimates + classification model" = {
                                           specific_vs_unspecific_table <- readxl::read_excel(file.path(task_file), sheet = "CM")
                                           setDT(specific_vs_unspecific_table)
                                           
                                           classification_rules_table <- readxl::read_excel(file.path(task_file), sheet = "Classification")
                                           setDT(classification_rules_table)
                                           
                                           endpoints_table <- readxl::read_excel(file.path(task_file), sheet = "Endpoints")
                                           setDT(endpoints_table)
                                           
                                           utils::write.table(specific_vs_unspecific_table, file = paste0(self$output_folder, '/SpecificVSUnspecific.csv'), sep=';', row.names = FALSE)
                                           utils::write.table(classification_rules_table, file = paste0(self$output_folder, '/ClassificationRules.csv'), sep=';', row.names = FALSE)                
                                           
                                           bmr_file <- readxl::read_excel(file.path(task_file), sheet = "BMRs")
                                           setDT(bmr_file)
                              
                                           self$fit_dt(experiments = 'merged', endpoints_table = endpoints_table)
                                           self$plotc(experiments = 'merged', for_endpoint = endpoints_table, bmr = bmr_file)
                                           self$classification_plot(experiments = 'merged', for_endpoint = endpoints_table, specific_vs_unspecific_table = specific_vs_unspecific_table, classification_rules_table = classification_rules_table)
                                         },
                                 )
                               }
                             }
                             
                             
                           )
)


'%push%' <- function(m, s) {
  s$protocol <- c(s$protocol, m)
  m
}


convert_table <- function(method = c('merged', 'single'), tabelle){
  # parameter
  method <- match.arg(method)
  
  result <- data.frame()
  tmp_tabelle <- tabelle
  setDT(tmp_tabelle)
  
  
  if(method == 'merged'){
    experimente <- unique(tmp_tabelle[, get('experiment')])
    doses <- sort(unique(tmp_tabelle[, get('dose')]))
    
    result <- data.frame('dose' = doses)
    
    N_experimente <- length(experimente)
    N_doses <- length(doses)
    for(e in 1:N_experimente){
      tmp_experiment <- experimente[e]
      
      tmp_col <- data.frame(rep(NA, N_doses))
      for(d in 1:N_doses){
        tmp_dose <- doses[d]
        val <- tmp_tabelle[(get('experiment') == tmp_experiment) & (get('dose') == tmp_dose), get('endpoint')]
        tmp_col[d, 1] <- ifelse(length(val)==0, as.numeric(NA), val)
      }
      colnames(tmp_col) <- paste0('experiment.', e)
      result <- cbind(result, tmp_col)
    }    
  }
  else if(method == 'single'){
    doses <- sort(unique(tmp_tabelle[, get('dose')]))
    
    N_doses <- length(doses)
    for(d in 1:N_doses){
      
      tmp_dose <- doses[d]
      vals <- t(tmp_tabelle[get('dose') == tmp_dose, get('endpoint')])
      
      tmp_row <- cbind(tmp_dose, vals)
      colnames(tmp_row)[1] <- 'dose'
      result <- rbind(result, tmp_row)
    }
  }
  
  result
}

# Returns string without leading or trailing white space
trim <- function (x) gsub("^\\s+|\\s+$", "", x)
