################################################################################
# Search database function -----------------------------------------------------
  # get_compound ---------------------------------------------------------------
#' @title get_compound
#' @author Zhiwei Zhou
#' @description search compound information via lab_id
#' @param lab_id Dodd lab ID
#' @importFrom magrittr %>%
#' @importFrom crayon blue red yellow green bgRed
#' @importFrom stringr str_detect str_extract
#' @export
#' @examples
#' \dontrun{
#' get_compound(lab_id = "S00001")
#' get_compound(lab_id = "MP000015")
#' get_compound(lab_id = "PP000015")
#' }


# get_compound(lab_id = "S00001")

setGeneric(name = "get_compound",
           def = function(lab_id){
             if (length(lab_id) > 1) {
               message(crayon::red("Please search one compound per time\n"))
               stop()
             }

             if (stringr::str_detect(lab_id, 'S\\d+|I\\d+')) {
               data('cpd_stanford_lib', envir = environment())
               data('cpd_iroa_lib', envir = environment())

               cpd_table <- cpd_stanford_lib$compound_table %>%
                 dplyr::bind_rows(cpd_iroa_lib$compound_table)
               rm(cpd_iroa_lib, cpd_stanford_lib);gc()

             }

             if (stringr::str_detect(lab_id, 'MP\\d+|MN\\d+')) {
               data('msdial_lib', envir = environment())

               cpd_table <- msdial_lib$positive$lib_meta %>%
                 dplyr::bind_rows(msdial_lib$negative$lib_meta) %>%
                 dplyr::rename(lab_id = id)

               rm(msdial_lib);gc()
             }

             if (stringr::str_detect(lab_id, 'PP\\d+|PN\\d+')) {
               data('fiehn_peptide_lib', envir = environment())

               cpd_table <- fiehn_peptide_lib$positive$lib_meta %>%
                 dplyr::bind_rows(fiehn_peptide_lib$negative$lib_meta) %>%
                 dplyr::rename(lab_id = id)

               rm(fiehn_peptide_lib);gc()
             }

             # match labid
             idx <- which(cpd_table$lab_id == lab_id)
             if (length(idx) == 0){
               message(crayon::red("This compound is not found in DoddLib\n"))
               stop()
             }
             if (length(idx) > 1) {
               message(crayon::red("This compound have multiple items in DoddLib\n"))
               stop()
             }

             cpd_table %>%
               dplyr::slice(idx) %>%
               dplyr::mutate_all(as.character) %>%
               tidyr::pivot_longer(everything()) %>%
               print(n = Inf)

           })


  # get_ms2 --------------------------------------------------------------------
#' @title get_ms2
#' @author Zhiwei Zhou
#' @description search compound ms2 spec via lab_id, ce, and polarity
#' @param lab_id Dodd lab ID
#' @param ce collision enery. "10", "20", "40"
#' @param polarity ionization polarity. "positive", "negative"
#' @importFrom magrittr %>%
#' @importFrom crayon blue red yellow green bgRed
#' @importFrom stringr str_detect str_extract
#' @export
#' @examples
#' \dontrun{
#' get_ms2(lab_id = 'S00001', ce = '10', polarity = 'positive')
#' get_ms2(lab_id = 'S00001', ce = '20', polarity = 'positive')
#' get_ms2(lab_id = 'I00101', ce = '20', polarity = 'negative')
#' get_ms2(lab_id = 'H00101', ce = '20', polarity = 'negative')
#' get_ms2(lab_id = 'MP000015', polarity = 'positive')
#' get_ms2(lab_id = 'PP000015', polarity = 'positive')
#' }

# get_ms2(lab_id = 'S00001', ce = '10', polarity = 'positive')
# get_ms2(lab_id = 'S00001', ce = '20', polarity = 'positive')
# get_ms2(lab_id = 'S00001', ce = '40', polarity = 'positive')
# get_ms2(lab_id = 'S00100', ce = '20', polarity = 'negative')
# get_ms2(lab_id = 'I00101', ce = '20', polarity = 'negative')
# get_ms2(lab_id = 'H00101', ce = '20', polarity = 'negative')
# get_ms2(lab_id = 'S00017', ce = '10', polarity = 'positive')
# get_ms2(lab_id = 'MP000015', polarity = 'positive')
# get_ms2(lab_id = 'PP000015', polarity = 'positive')

setGeneric(name = 'get_ms2',
           def = function(lab_id,
                          ce = c('10', '20', '40'),
                          polarity = c('positive', "negative")){

             ce <- match.arg(ce)
             polarity <- match.arg(polarity)


             # dodd lib
             if (stringr::str_detect(lab_id, 'S\\d+|I\\d+')) {
               data("ms2_dodd_lib", envir = environment())
               data("ms2_stanford_lib", envir = environment())
               data("ms2_iroa_lib", envir = environment())

               ids_pos <- names(ms2_dodd_lib$positive)
               ids_neg <- names(ms2_dodd_lib$positive)

               # search dodd lib first, if not, according id to search separated db
               if (polarity == 'positive' & (lab_id %in% ids_pos)) {
                 temp_spec <- ms2_dodd_lib$positive[[lab_id]][[ce]]
                 if (is.null(temp_spec)) {
                   message(crayon::red("No spectra of this compound & this ce is found!"))
                 } else {
                   return(temp_spec)
                 }
               }

               if (polarity == 'negative' & (lab_id %in% ids_neg)) {
                 temp_spec <- ms2_dodd_lib$negative[[lab_id]][[ce]]
                 if (is.null(temp_spec)) {
                   message(crayon::red("No spectra of this compound & this ce is found!"))
                 } else {
                   return(temp_spec)
                 }
               }

               if (stringr::str_detect(lab_id, pattern = "S\\d+")) {
                 if (polarity == "positive") {
                   temp_spec1 <- ms2_stanford_lib$hilic_pos[[lab_id]][[ce]]
                   temp_spec2 <- ms2_stanford_lib$c18_pos[[lab_id]][[ce]]

                   # if two spec found, return the highest one
                   if (all(!is.null(temp_spec1), !is.null(temp_spec2))) {

                     if (nrow(temp_spec1) >= 10) {
                       temp_int1 <- temp_spec1[,2] %>% sort(decreasing = TRUE) %>% .[1:10] %>% sum()
                     } else {
                       temp_int1 <- temp_spec1[,2] %>% sum()
                     }

                     if (nrow(temp_spec2) >= 10) {
                       temp_int2 <- temp_spec2[,2] %>% sort(decreasing = TRUE) %>% .[1:10] %>% sum()
                     } else {
                       temp_int2 <- temp_spec2[,2] %>% sum()
                     }

                     if (temp_int1 > temp_int2) {
                       message(crayon::yellow("Two spectra of this compound are found, the highest spectrum (HILIC_Pos) is returned!"))
                       return(temp_spec1)
                     } else {
                       message(crayon::yellow("Two spectra of this compound are found, the highest spectrum (C18_Pos) is returned!"))
                       return(temp_spec2)
                     }

                   }

                   if (!is.null(temp_spec1) & is.null(temp_spec2)) {
                     return(temp_spec1)
                   }

                   if (is.null(temp_spec1) & !is.null(temp_spec2)) {
                     return(temp_spec2)
                   }

                   if (is.null(temp_spec1) & is.null(temp_spec2)) {
                     message(crayon::red("No spectra of this compound is found!"))
                   }

                 } else {
                   temp_spec1 <- ms2_stanford_lib$hilic_neg[[lab_id]][[ce]]
                   temp_spec2 <- ms2_stanford_lib$c18_neg[[lab_id]][[ce]]

                   # if two spec found, return the highest one
                   if (all(!is.null(temp_spec1), !is.null(temp_spec2))) {

                     if (nrow(temp_spec1) >= 10) {
                       temp_int1 <- temp_spec1[,2] %>% sort(decreasing = TRUE) %>% .[1:10] %>% sum()
                     } else {
                       temp_int1 <- temp_spec1[,2] %>% sum()
                     }

                     if (nrow(temp_spec2) >= 10) {
                       temp_int2 <- temp_spec2[,2] %>% sort(decreasing = TRUE) %>% .[1:10] %>% sum()
                     } else {
                       temp_int2 <- temp_spec2[,2] %>% sum()
                     }

                     if (temp_int1 > temp_int2) {
                       message(crayon::yellow("Two spectra of this compound are found, the highest spectrum (HILIC_Pos) is returned!"))
                       return(temp_spec1)
                     } else {
                       message(crayon::yellow("Two spectra of this compound are found, the highest spectrum (C18_Pos) is returned!"))
                       return(temp_spec2)
                     }

                   }

                   if (!is.null(temp_spec1) & is.null(temp_spec2)) {
                     return(temp_spec1)
                   }

                   if (is.null(temp_spec1) & !is.null(temp_spec2)) {
                     return(temp_spec2)
                   }

                   if (is.null(temp_spec1) & is.null(temp_spec2)) {
                     message(crayon::red("No spectra of this compound is found!"))
                   }
                 }
               } else if (stringr::str_detect(lab_id, pattern = "I\\d+")) {
                 if (polarity == "positive") {
                   temp_spec1 <- ms2_iroa_lib$hilic_pos[[lab_id]][[ce]]
                   temp_spec2 <- ms2_iroa_lib$c18_pos[[lab_id]][[ce]]

                   # if two spec found, return the highest one
                   if (all(!is.null(temp_spec1), !is.null(temp_spec2))) {

                     if (nrow(temp_spec1) >= 10) {
                       temp_int1 <- temp_spec1[,2] %>% sort(decreasing = TRUE) %>% .[1:10] %>% sum()
                     } else {
                       temp_int1 <- temp_spec1[,2] %>% sum()
                     }

                     if (nrow(temp_spec2) >= 10) {
                       temp_int2 <- temp_spec2[,2] %>% sort(decreasing = TRUE) %>% .[1:10] %>% sum()
                     } else {
                       temp_int2 <- temp_spec2[,2] %>% sum()
                     }

                     if (temp_int1 > temp_int2) {
                       message(crayon::yellow("Two spectra of this compound are found, the highest spectrum (HILIC_Pos) is returned!"))
                       return(temp_spec1)
                     } else {
                       message(crayon::yellow("Two spectra of this compound are found, the highest spectrum (C18_Pos) is returned!"))
                       return(temp_spec2)
                     }

                   }

                   if (!is.null(temp_spec1) & is.null(temp_spec2)) {
                     return(temp_spec1)
                   }

                   if (is.null(temp_spec1) & !is.null(temp_spec2)) {
                     return(temp_spec2)
                   }

                   if (is.null(temp_spec1) & is.null(temp_spec2)) {
                     message(crayon::red("No spectra of this compound is found!"))
                   }

                 } else {
                   temp_spec1 <- ms2_iroa_lib$hilic_neg[[lab_id]][[ce]]
                   temp_spec2 <- ms2_iroa_lib$c18_neg[[lab_id]][[ce]]

                   # if two spec found, return the highest one
                   if (all(!is.null(temp_spec1), !is.null(temp_spec2))) {

                     if (nrow(temp_spec1) >= 10) {
                       temp_int1 <- temp_spec1[,2] %>% sort(decreasing = TRUE) %>% .[1:10] %>% sum()
                     } else {
                       temp_int1 <- temp_spec1[,2] %>% sum()
                     }

                     if (nrow(temp_spec2) >= 10) {
                       temp_int2 <- temp_spec2[,2] %>% sort(decreasing = TRUE) %>% .[1:10] %>% sum()
                     } else {
                       temp_int2 <- temp_spec2[,2] %>% sum()
                     }

                     if (temp_int1 > temp_int2) {
                       message(crayon::yellow("Two spectra of this compound are found, the highest spectrum (HILIC_Pos) is returned!"))
                       return(temp_spec1)
                     } else {
                       message(crayon::yellow("Two spectra of this compound are found, the highest spectrum (C18_Pos) is returned!"))
                       return(temp_spec2)
                     }

                   }

                   if (!is.null(temp_spec1) & is.null(temp_spec2)) {
                     return(temp_spec1)
                   }

                   if (is.null(temp_spec1) & !is.null(temp_spec2)) {
                     return(temp_spec2)
                   }

                   if (is.null(temp_spec1) & is.null(temp_spec2)) {
                     message(crayon::red("No spectra of this compound is found!"))
                   }
                 }
               } else {
                 message(crayon::red("No spectra of this compound is found!"))
               }
             }


             if (stringr::str_detect(lab_id, 'MP\\d+|MN\\d+')) {
               data("msdial_lib", envir = environment())
               temp_spec <- msdial_lib[[polarity]][['lib_spec']][[lab_id]]
               rm(msdial_lib);gc()
               if (is.null(temp_spec)) {
                 message(crayon::red("No spectra of this compound & this ce is found!"))
               } else {
                 return(temp_spec)
               }
             }

             if (stringr::str_detect(lab_id, 'PP\\d+|PN\\d+')) {
               data("fiehn_peptide_lib", envir = environment())
               temp_spec <- fiehn_peptide_lib[[polarity]][['lib_spec']][[lab_id]]
               rm(fiehn_peptide_lib);gc()
               if (is.null(temp_spec)) {
                 message(crayon::red("No spectra of this compound & this ce is found!"))
               } else {
                 return(temp_spec)
               }
             }

           })


