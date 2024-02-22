################################################################################
# load_spec_db -----------------------------------------------------------------

#' @title load_spec_db
#' @author Zhiwei Zhou
#' @param lib database name, 'dodd', 'msdial', 'peptide'. Default: 'dodd'
#' @param column 'hilic', 'c18'. Default: 'hilic'
#' @param ce "10", "20", "40"; Default: '20'
#' @param polarity 'positive' or 'negative'. Default: 'positive'
#' @param adduct_list NULL
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
#' # peptide lib
#' test <- load_spec_db(lib = 'peptide', polarity = 'positive')
#' }
#'

# test <- load_spec_db(lib = 'dodd', column = 'hilic', ce = '20', polarity = 'positive', adduct_list = '[M+H]+')
# test <- load_spec_db(lib = 'msdial', polarity = 'positive')
# test <- load_spec_db(lib = 'peptide', polarity = 'positive')
# test <- load_spec_db(lib = 'gnps_bile_acid', polarity = 'positive')
# test <- load_spec_db(lib = 'gnps_acyl_amides', polarity = 'positive')
# test <- load_spec_db(lib = 'gnps_acyl_esters', polarity = 'positive')
# test <- load_spec_db(lib = 'all_public', polarity = 'positive')

setGeneric(name = 'load_spec_db',
           def = function(
    lib = c('dodd', 'msdial', 'peptide', 'gnps_bile_acid', 'gnps_acyl_amides', 'gnps_acyl_esters', 'all_public'),
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
                   },
                   'msdial' = {
                     paste0('Load the ', lib, ' library\n',
                            'Parameters:\n',
                            'Polarity: ', polarity, '\n')
                   },
                   'gnps_bile_acid' = {
                     paste0('Load the ', lib, ' library\n',
                            'Parameters:\n',
                            'Polarity: ', polarity, '\n')
                   },
                   'gnps_acyl_amides' = {
                     paste0('Load the ', lib, ' library\n',
                            'Parameters:\n',
                            'Polarity: ', polarity, '\n')
                   },
                   'gnps_acyl_esters' = {
                     paste0('Load the ', lib, ' library\n',
                            'Parameters:\n',
                            'Polarity: ', polarity, '\n')
                   },
                   'all_public' = {
                     paste0('Load the ', lib, ' library\n',
                            'Parameters:\n',
                            'Polarity: ', polarity, '\n')
                   },
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
                     },
                     'peptide' = {
                       data("fiehn_peptide_lib", envir = environment())
                       cpd_lib <- fiehn_peptide_lib[[polarity]]$lib_meta
                       ms2_lib <- fiehn_peptide_lib[[polarity]]$lib_spec
                       rm(fiehn_peptide_lib);gc()
                     },
                     'gnps_bile_acid' = {
                       data("gnps_bile_acid_lib", envir = environment())
                       cpd_lib <- gnps_bile_acid_lib[[polarity]]$lib_meta
                       ms2_lib <- gnps_bile_acid_lib[[polarity]]$lib_spec
                       rm(gnps_bile_acid_lib);gc()
                     },
                     'gnps_acyl_amides' = {
                       data("gnps_acyl_amides_lib", envir = environment())
                       cpd_lib <- gnps_acyl_amides_lib[[polarity]]$lib_meta
                       ms2_lib <- gnps_acyl_amides_lib[[polarity]]$lib_spec
                       rm(gnps_acyl_amides_lib);gc()
                     },
                     'gnps_acyl_esters' = {
                       data("gnps_acyl_esters_lib", envir = environment())
                       cpd_lib <- gnps_acyl_esters_lib[[polarity]]$lib_meta
                       ms2_lib <- gnps_acyl_esters_lib[[polarity]]$lib_spec
                       rm(gnps_acyl_esters_lib);gc()
                     },
                     'all_public' = {
                       data("msdial_lib", envir = environment())
                       data("fiehn_peptide_lib", envir = environment())
                       data("gnps_bile_acid_lib", envir = environment())
                       data("gnps_acyl_amides_lib", envir = environment())
                       data("gnps_acyl_esters_lib", envir = environment())

                       external_id_msdial <- msdial_lib[[polarity]]$lib_meta$comment %>%
                         stringr::str_extract('DB#=.+;') %>%
                         stringr::str_replace('DB#=', '') %>%
                         stringr::str_replace(';', '')

                       external_id_bile_acid <- gnps_bile_acid_lib[[polarity]]$lib_meta$comment %>%
                         stringr::str_extract('DB#=.+;') %>%
                         stringr::str_replace('DB#=', '') %>%
                         stringr::str_replace(';', '')

                       external_id_acyl_amides <- gnps_acyl_amides_lib[[polarity]]$lib_meta$comment %>%
                         stringr::str_extract('DB#=.+;') %>%
                         stringr::str_replace('DB#=', '') %>%
                         stringr::str_replace(';', '')

                       external_id_acyl_esters <- gnps_acyl_esters_lib[[polarity]]$lib_meta$comment %>%
                         stringr::str_extract('DB#=.+;') %>%
                         stringr::str_replace('DB#=', '') %>%
                         stringr::str_replace(';', '')

                       idx_rep_bile_acid <- which(external_id_bile_acid %in% external_id_msdial)
                       idx_rep_acyl_amides <- which(external_id_acyl_amides %in% external_id_msdial)
                       idx_rep_acyl_esters <- which(external_id_acyl_esters %in% external_id_msdial)

                       cpd_lib1 <- msdial_lib[[polarity]]$lib_meta %>% dplyr::select(id:inchikey)
                       ms2_lib1 <- msdial_lib[[polarity]]$lib_spec
                       cpd_lib2 <- fiehn_peptide_lib[[polarity]]$lib_meta  %>% dplyr::select(id:inchikey)
                       ms2_lib2 <- fiehn_peptide_lib[[polarity]]$lib_spec

                       cpd_lib3 <- gnps_bile_acid_lib[[polarity]]$lib_meta %>% dplyr::select(id:inchikey)
                       ms2_lib3 <- gnps_bile_acid_lib[[polarity]]$lib_spec
                       if (length(idx_rep_bile_acid) > 0) {
                         cpd_lib3 <- cpd_lib3[-idx_rep_bile_acid,]
                         ms2_lib3 <- ms2_lib3[-idx_rep_bile_acid]
                       }
                       cpd_lib4 <- gnps_acyl_amides_lib[[polarity]]$lib_meta %>% dplyr::select(id:inchikey)
                       ms2_lib4 <- gnps_acyl_amides_lib[[polarity]]$lib_spec
                       if (length(idx_rep_acyl_amides) > 0) {
                         cpd_lib4 <- cpd_lib4[-idx_rep_acyl_amides,]
                         ms2_lib4 <- ms2_lib4[-idx_rep_acyl_amides]
                       }
                       cpd_lib5 <- gnps_acyl_esters_lib[[polarity]]$lib_meta %>% dplyr::select(id:inchikey)
                       ms2_lib5 <- gnps_acyl_esters_lib[[polarity]]$lib_spec
                       if (length(idx_rep_acyl_esters) > 0) {
                         cpd_lib5 <- cpd_lib4[-idx_rep_acyl_esters,]
                         ms2_lib5 <- ms2_lib4[-idx_rep_acyl_esters]
                       }

                       cpd_lib <- cpd_lib1 %>%
                         dplyr::bind_rows(cpd_lib2) %>%
                         dplyr::bind_rows(cpd_lib3) %>%
                         dplyr::bind_rows(cpd_lib4) %>%
                         dplyr::bind_rows(cpd_lib5)

                       ms2_lib <- c(ms2_lib1,
                                    ms2_lib2,
                                    ms2_lib3,
                                    ms2_lib4,
                                    ms2_lib5)

                       rm(msdial_lib, fiehn_peptide_lib, gnps_bile_acid_lib, gnps_acyl_amides_lib, gnps_acyl_esters_lib);gc()
                       rm(cpd_lib1, cpd_lib2, cpd_lib3, cpd_lib4, cpd_lib5);gc()
                       rm(ms2_lib1, ms2_lib2, ms2_lib3, ms2_lib4, ms2_lib5);gc()
                       rm(external_id_msdial, external_id_bile_acid, external_id_acyl_amides, external_id_acyl_esters);gc()
                       rm(idx_rep_bile_acid, idx_rep_acyl_amides, idx_rep_acyl_esters);gc()
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
               if (length(adduct_list) > 0) {
                 cpd_lib <- cpd_lib %>% dplyr::filter(adduct %in% adduct_list)
                 if (nrow(cpd_lib) == 0) {
                   stop('Please modify the adduct form because no record is available')
                 }
                 lib_spec <- match(cpd_lib$id, names(ms2_lib)) %>% ms2_lib[.]
               }

               result <- list(lib_meta = cpd_lib,
                              lib_spec = ms2_lib)

               return(result)
             }

             # fiehn peptide lib -----------------------------------------------
             if (lib == 'peptide') {
               if (length(adduct_list) > 0) {
                 cpd_lib <- cpd_lib %>% dplyr::filter(adduct %in% adduct_list)
                 if (nrow(cpd_lib) == 0) {
                   stop('Please modify the adduct form because no record is available')
                 }
                 lib_spec <- match(cpd_lib$id, names(ms2_lib)) %>% ms2_lib[.]
               }

               result <- list(lib_meta = cpd_lib,
                              lib_spec = ms2_lib)

               return(result)
             }

             # gnps_bile_acid lib -----------------------------------------------
             if (lib == 'gnps_bile_acid') {
               if (length(adduct_list) > 0) {
                 cpd_lib <- cpd_lib %>% dplyr::filter(adduct %in% adduct_list)
                 if (nrow(cpd_lib) == 0) {
                   stop('Please modify the adduct form because no record is available')
                 }
                 lib_spec <- match(cpd_lib$id, names(ms2_lib)) %>% ms2_lib[.]
               }
               result <- list(lib_meta = cpd_lib,
                              lib_spec = ms2_lib)

               return(result)
             }

             # gnps_acyl_amides lib -----------------------------------------------
             if (lib == 'gnps_acyl_amides') {
               if (length(adduct_list) > 0) {
                 cpd_lib <- cpd_lib %>% dplyr::filter(adduct %in% adduct_list)
                 if (nrow(cpd_lib) == 0) {
                   stop('Please modify the adduct form because no record is available')
                 }
                 lib_spec <- match(cpd_lib$id, names(ms2_lib)) %>% ms2_lib[.]
               }
               result <- list(lib_meta = cpd_lib,
                              lib_spec = ms2_lib)

               return(result)
             }

             # gnps_acyl_amides lib -----------------------------------------------
             if (lib == 'gnps_acyl_esters') {
               if (length(adduct_list) > 0) {
                 cpd_lib <- cpd_lib %>% dplyr::filter(adduct %in% adduct_list)
                 if (nrow(cpd_lib) == 0) {
                   stop('Please modify the adduct form because no record is available')
                 }
                 lib_spec <- match(cpd_lib$id, names(ms2_lib)) %>% ms2_lib[.]
               }
               result <- list(lib_meta = cpd_lib,
                              lib_spec = ms2_lib)

               return(result)
             }

             # gnps_acyl_amides lib -----------------------------------------------
             if (lib == 'gnps_acyl_esters') {
               if (length(adduct_list) > 0) {
                 cpd_lib <- cpd_lib %>% dplyr::filter(adduct %in% adduct_list)
                 if (nrow(cpd_lib) == 0) {
                   stop('Please modify the adduct form because no record is available')
                 }
                 lib_spec <- match(cpd_lib$id, names(ms2_lib)) %>% ms2_lib[.]
               }
               result <- list(lib_meta = cpd_lib,
                              lib_spec = ms2_lib)

               return(result)
             }

             # gnps_acyl_amides lib -----------------------------------------------
             if (lib == 'all_public') {
               if (length(adduct_list) > 0) {
                 cpd_lib <- cpd_lib %>% dplyr::filter(adduct %in% adduct_list)
                 if (nrow(cpd_lib) == 0) {
                   stop('Please modify the adduct form because no record is available')
                 }
                 lib_spec <- match(cpd_lib$id, names(ms2_lib)) %>% ms2_lib[.]
               }
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



################################################################################
# startup massage --------------------------------------------------------------
.onAttach <- function(libname, pkgname){
  packageStartupMessage("
Version 0.2.2
-------------
Authors: Zhiwei Zhou
Maintainer: Zhiwei Zhou

Updates
-------------
o Add pathdb_enteropathway database
o Add startup message
")
}
