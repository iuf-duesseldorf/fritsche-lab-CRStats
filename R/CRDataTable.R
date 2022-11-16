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
                               N <- nrow(self$dt)
                               for(i in 1:N){
                                 
                                 ccounter <- self$dt[i, get(counter_col)]
                                 cdenominator <- self$dt[i, get(denominator_col)]
                                 
                                 if(is.na(ccounter) || is.na(cdenominator)){
                                   
                                   self$dt[i, new_col_name_4_relative] <- as.numeric(NA)
                                   
                                 }
                                 else if(cdenominator == 0){
                                   
                                   self$dt[i, new_col_name_4_relative] <- 0
                                   
                                 }
                                 
                               }
                               
                               
                               # is.na(self$dt) <- sapply(self$dt, is.infinite)
                               # self$dt[, (new_col_name_4_relative) := lapply(.SD, function(x) ifelse(is.infinite(x), 0, x)), .SD = new_col_name_4_relative] # x/0 case
                               # self$dt[, (new_col_name_4_relative) := lapply(.SD, function(x) ifelse(is.nan(x), 0, x)), .SD = new_col_name_4_relative] # 0/0 case
                               
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
                                 self$save_dt('AddedQuotients')
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
                                   
                                   tmp_dose <- drc_model[["model"]][["dose"]]
                                   min_index <- which.min(tmp_dose) # we do not want to require ordered dataset
                                   
                                   if(best_fit$model == 'no-effect') min_index <- 1;
                                   
                                   lowest_dose_predres <- drc_model[["fitted.values"]][[min_index]]
                                   
                                 }
                                 else{
                                   
                                   tmp_dose <- drc_model[["data"]][["dose"]]
                                   min_index <- which.min(tmp_dose) # we do not want to require ordered dataset
                                   lowest_dose_predres <- drc_model[["predres"]][, 'Predicted values'][[min_index]]
                                   
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
                             }
                             
                           )
)


'%push%' <- function(m, s) {
  s$protocol <- c(s$protocol, m)
  m
}
