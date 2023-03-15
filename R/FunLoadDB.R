################################################################################
# load_spec_db -----------------------------------------------------------------

#' @title load_spec_db
#' @author Zhiwei Zhou
#' @param lib database name, 'zhuMetLib', 'zhuMetLib_orbitrap', 'fiehnHilicLib'. Default: 'zhuMetLib'
#' @param instrument "SciexTripleTOF", "AgilentQTOF", "BrukerQTOF", "ThermoOrbitrap", "ThermoExploris", "AgilentDTIMMS", "BrukerTIMS". Default: "SciexTripleTOF"
#' @param column 'hilic', 'rp'. Default: 'hilic'
#' @param method_lc 'Amide12min' or 'Amide23min'
#' @param ce "10", "20", "30", "35,15", "40", "50"; Default: '30'
#' @param polarity 'positive' or 'negative'. Default: 'positive'
#' @param adduct_list NULL
# #' @param is_rt_score whether only reserve compounds with RT
# #' @param is_ccs_score whether only reserve compounds with CCS
#' @param is_rt_calibration TRUE
#' @param path '.'
#' @importFrom magrittr %>%
#' @importFrom crayon blue red yellow green bgRed
#' @importFrom stringr str_detect str_extract
#' @export
#' @examples
#' \dontrun{
#' # dodd lib
#' test <- load_spec_db(lib = 'dodd', column = 'hilic', ce = '20', polarity = 'positive', adduct_list = '[M+H]+')
#' # msdial lib
#' test <- load_spec_db(lib = 'msdial', polarity = 'positive')
#' }
#'

# test <- load_spec_db(lib = 'dodd', column = 'hilic', ce = '20', polarity = 'positive', adduct_list = '[M+H]+')
# test <- load_spec_db(lib = 'msdial', polarity = 'positive')

setGeneric(name = 'load_spec_db',
           def = function(
    lib = c('dodd', 'msdial'),
    column = c('hilic', 'c18'),
    ce = c('10', '20', '40'),
    polarity = c('positive', 'negative'),
    adduct_list = NULL
    # file_rt_ref = 'RT_recalibration_table.csv',
    # is_rt_calibration = FALSE,
    # path = NULL
           ){
             lib <- match.arg(lib)
             polarity <- match.arg(polarity)

             message(
               crayon::blue(
                 switch (lib,
                   'dodd' = {
                     paste0('Load the ', lib, ' library\n',
                            'Parameters:\n',
                            'Collision energy: ', ce, '\n',
                            'Column: ', column, '\n',
                            'Polarity: ', polarity, '\n')
                   },
                   'msdial' = {
                     paste0('Load the ', lib, ' library\n',
                            'Parameters:\n',
                            'Polarity: ', polarity, '\n')
                   }
                 )
               )
             )

             switch (lib,
                     'dodd' = {
                       data('cpd_dodd_lib', envir = environment())
                       data('ms2_dodd_lib', envir = environment())
                       cpd_lib <- cpd_dodd_lib
                       ms2_lib <- ms2_dodd_lib
                       rm(cpd_dodd_lib, ms2_dodd_lib);gc()
                     },
                     'msdial' = {
                       data("msdial_lib", envir = environment())
                       cpd_lib <- msdial_lib[[polarity]]$lib_meta
                       ms2_lib <- msdial_lib[[polarity]]$lib_spec
                       rm(msdial_lib);gc()
                     }
             )

             # dodd lib --------------------------------------------------------
             if (lib == 'dodd') {
               column <- match.arg(column)
               ce <- match.arg(ce)
               # check adduct type
               if (length(adduct_list) == 0) {
                 stop('Please input adduct_list\n')
               }

               if (polarity == 'positive') {
                 temp <- all(adduct_list %in% lib_adduct_nl$positive$adduct)
                 if (!temp) stop('Please check the selected adducts!\n')
               } else {
                 temp <- all(adduct_list %in% lib_adduct_nl$negative$adduct)
                 if(!temp) stop('Please check the selected adducts!\n')
               }



               # modify lib_meta
               lib_meta <- cpd_lib$compound_table
               switch(polarity,
                      'positive' = {
                        temp_adduct_table <- lib_adduct_nl$positive %>% dplyr::filter(adduct %in% adduct_list)
                      },
                      'negative' = {
                        temp_adduct_table <- lib_adduct_nl$negative %>% dplyr::filter(adduct %in% adduct_list)
                      }
               )

               # calculate mz
               lib_mz <- lapply(seq_along(lib_meta$lab_id), function(i){
                 temp_mz <- lib_meta$monoisotopic_mass[i]

                 result <- sapply(seq_along(temp_adduct_table$adduct), function(j){
                   calculateMz(exact_mass = temp_mz,
                               adduct = temp_adduct_table$adduct[j],
                               delta_mz = temp_adduct_table$delta_mz[j])
                 })
                 result
               })
               lib_mz <- lib_mz %>% do.call(rbind, .) %>% tibble::as_tibble()
               colnames(lib_mz) <- temp_adduct_table$adduct

               # extract compound info
               lib_meta <- lib_meta %>%
                 dplyr::bind_cols(lib_mz) %>%
                 tidyr::pivot_longer(cols = temp_adduct_table$adduct,
                                     names_to = 'adduct',
                                     values_to = 'mz') %>%
                 dplyr::rename(id = lab_id, name = compound_name) %>%
                 dplyr::select(id:formula, adduct:mz, dplyr::everything())

               # select mode responding RT
               if (column == 'hilic') {
                 if (polarity == 'positive') {
                   lib_meta <- lib_meta %>%
                     dplyr::mutate(rt = rt_hilic_pos) %>%
                     dplyr::mutate(rt = round(rt*60)) %>%
                     dplyr::select(id:mz, rt, dplyr::everything())
                 } else {
                   lib_meta <- lib_meta %>%
                     dplyr::mutate(rt = rt_hilic_neg) %>%
                     dplyr::mutate(rt = round(rt*60)) %>%
                     dplyr::select(id:mz, rt, dplyr::everything())
                 }
               }

               if (column == 'c18') {
                 if (polarity == 'positive') {
                   lib_meta <- lib_meta %>%
                     dplyr::mutate(rt = rt_c18_pos) %>%
                     dplyr::mutate(rt = round(rt*60)) %>%
                     dplyr::select(id:mz, rt, dplyr::everything())
                 } else {
                   lib_meta <- lib_meta %>%
                     dplyr::mutate(rt = rt_c18_neg) %>%
                     dplyr::mutate(rt = round(rt*60)) %>%
                     dplyr::select(id:mz, rt, dplyr::everything())
                 }
               }

               # modify lib_spec
               if (polarity == 'positive') {
                 lib_spec <- ms2_lib$positive
               } else {
                 lib_spec <- ms2_lib$negative
               }

               idx_spec <- match(unique(lib_meta$id), names(lib_spec))
               idx_spec <- which(!is.na(idx_spec))
               ms2_id <- unique(lib_meta$id)[idx_spec]
               lib_spec <- match(ms2_id, names(lib_spec)) %>%
                 lib_spec[.] %>%
                 lapply(., function(x){
                   x[[ce]]
                 })

               # rm null spec
               idx_null <- sapply(lib_spec, is.null) %>% which()
               lib_spec <- lib_spec[-idx_null]
               result <- list(lib_meta = lib_meta,
                              lib_spec = lib_spec)

               return(result)
             }

             # msdial lib ------------------------------------------------------
             if (lib == 'msdial') {
               result <- list(lib_meta = cpd_lib,
                              lib_spec = ms2_lib)

               return(result)
             }



           }
)




#   calibrateRT ----------------------------------------------------------------

#' @title calibrateRT
#' @description RT calibration according to the RTQC
#' @author Zhiwei Zhou, Mingdu Luo
#' \email{zhouzw@@sioc.ac.cn}
#' @param file_rt_ref Default: 'RT_recalibration_table.csv'
#' @param lib_rt a vector of RT values
#' @param is_rt_calibration whether need to calibrate RT
#' @param is_plot Default: TRUE
#' @param method_lc 'Amide12min', 'Amide23min'. Default:  'Amide12min'
#' @param column Default: 'hilic'
#' @param path '.'

setGeneric(name = 'calibrateRT',
           def = function(file_rt_ref = 'RT_recalibration_table.csv',
                          lib_rt,
                          is_rt_calibration = TRUE,
                          is_plot = TRUE,
                          method_lc = c('Amide12min', 'Amide23min'),
                          column = 'hilic',
                          path = '.'){

             # Load in house RT library according to LC method
             switch(method_lc,
                    "Amide12min" = {
                      data('rt_ref', envir = environment())
                      ref_rtqc_table <- rt_ref[[1]]
                      lc_start <- 0
                      lc_end <- 720},

                    "Amide23min" = {
                      data('rt_ref', envir = environment())
                      ref_rtqc_table <- rt_ref[[2]]
                      lc_start <- 0
                      lc_end <- 1380}
             )

             # Check whether to do RT recalibration, only hilic method were allowed ###
             if (any(!is_rt_calibration, column != 'hilic')){
               cat('RT recalibration was turned off.\n')
             } else {

               ### Check existence of RT calibration table ###
               if (!("RT_recalibration_table.csv" %in% dir(path))){
                 stop("There was no 'RT_recalibration_table.csv' table in the file folder.\n")
               } else {

                 exp_rtqc_rt <- read.csv(file.path(path, file_rt_ref), stringsAsFactors = F)

                 ### Check the format of exp_rtqc_rt ###
                 if (!identical(colnames(exp_rtqc_rt)[1:6], c("compound.name", "id.zhulab", "id.pubchem", "ref.mz", "rt", "polarity"))){
                   stop("Please check the format of 'RT_recalibration_table.csv' table and correct it according to our tutorial.\n")
                 } else {

                   ### Begin RT calibration ###
                   if (max(exp_rtqc_rt$rt) < 60) {
                     exp_rtqc_rt$rt <- round(exp_rtqc_rt$rt*60, digits = 4)
                   } else {
                     exp_rtqc_rt$rt <- round(exp_rtqc_rt$rt, digits = 4)
                   }

                   # the index of match rt compound and the warning
                   idx <- match(toupper(ref_rtqc_table$name),
                                toupper(exp_rtqc_rt$id.zhulab))

                   if (sum(!is.na(idx)) < 7){
                     warning(paste0("The number of used compouds for RT recalibration with LOESS was ",
                                    sum(!is.na(idx)),
                                    " and might be insufficient for a good performance.\n\n"))
                   } else {
                     cat("The number of used compouds for RT recalibration with LOESS was ",
                         sum(!is.na(idx)),
                         ".\n\n",
                         sep ='')
                   }

                   training.data <- data.frame(ref.rt = c(lc_start, ref_rtqc_table$rt[!is.na(idx)], lc_end),
                                               exp.rt = c(lc_start, exp_rtqc_rt$rt[idx[!is.na(idx)]], lc_end))

                   rownames(training.data) <- c('Start', exp_rtqc_rt$id.zhulab[idx[!is.na(idx)]], 'End')

                   ### rt.recalibration.model <- lm(exp.rt~ref.rt, data = training.data)
                   rt.recalibration.model <- loess(exp.rt~ref.rt,
                                                   data = training.data,
                                                   span = 0.75, degree = 2)

                   new.data <- data.frame(ref.rt = lib_rt$rt, stringsAsFactors = FALSE)
                   lib.calibrated.rt <- round(predict(object =  rt.recalibration.model,
                                                      newdata = new.data),
                                              digits = 2)

                   lib_rt$rt <- lib.calibrated.rt

                   result <- list(lib_rt = lib_rt,
                                  training.data = training.data,
                                  rt.recalibration.model = rt.recalibration.model)

                   dir.create(file.path(path, "01_result_initial_seed_annotation", '00_intermediate_data'), showWarnings = FALSE, recursive = TRUE)
                   save(result,
                        file = file.path(path, "01_result_initial_seed_annotation", '00_intermediate_data', 'rt_calibration_result'),
                        version = 2)

                   cat("RT recalibration was done.\n")

                   if (is_plot) {
                     # browser()
                     dir.create(file.path(path, "01_result_initial_seed_annotation", '01_rt_calibration_plot'),
                                showWarnings = FALSE, recursive = TRUE)
                     plotRtCalibration(result_rt_calibration = result[[2]],
                                       rt_recalibration_model = result[[3]],
                                       path = file.path(path, "01_result_initial_seed_annotation", '01_rt_calibration_plot'))
                   }

                   return(result)
                   rm(c(idx, lc_start, lc_end))

                 }
               }
             }
           })





#   calculateMz ----------------------------------------------------------------
setGeneric(name = 'calculateMz',
           def = function(
    exact_mass,
    adduct,
    delta_mz,
    nmol = NULL,
    ncharge = NULL
           ){

             if (length(nmol) == 0) {
               if (stringr::str_detect(adduct, pattern = '2M')) {
                 mz <- exact_mass*2 + delta_mz
               } else if (stringr::str_detect(adduct, pattern = '3M')) {
                 mz <- exact_mass*3 + delta_mz
               } else {
                 mz <- exact_mass + delta_mz
               }
             } else {
               mz <- exact_mass*nmol + delta_mz
             }


             if (length(ncharge) == 0) {
               if (stringr::str_detect(adduct, pattern = '\\]2\\-|\\]2\\+')) {
                 mz <- mz/2
               } else if (stringr::str_detect(adduct, pattern = '\\]3\\-|\\]3\\+')) {
                 mz <- mz/3
               } else {
                 mz
               }
             } else {
               mz <- mz/ncharge
             }

             mz
           })
