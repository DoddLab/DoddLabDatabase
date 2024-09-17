#' @title export_db_tidymass
#' @author Zhiwei Zhou
#' @param lib dodd lab library. "stanford", "iroa"
#' @param column columne type. "hilic", "rp"
#' @param polarity ionization polarity. "positive", "negative"
#' @importFrom magrittr %>%
#' @importFrom crayon blue red yellow green bgRed
#' @importFrom stringr str_detect str_extract
#' @export
#' @examples
#' \dontrun{
#' stanford_hilic_pos_v1.0 <- export_db_tidymass(lib = 'stanford',
#' column = 'hilic',
#' polarity = 'positive')
#' }

# stanford_hilic_pos_v1.0 <- export_db_tidymass(lib = 'stanford',
#                                               column = 'hilic',
#                                               polarity = 'positive')


# Dodd lab has different gradients for positive/negative modes, so the database should be exported specific to column, polarity

setGeneric(name = 'export_db_tidymass',
           def = function(lib = c('dodd', 'stanford', 'iroa'),
                          column = c('hilic', 'rp'),
                          polarity = c('positive', 'negative')){

             message(crayon::blue("Start export database for tidymass..."))

             lib <- match.arg(lib)
             column <- match.arg(column)
             polarity <- match.arg(polarity)

             # load databases
             switch (lib,
                     "dodd" = {
                       data('cpd_dodd_lib', envir = environment())
                       data('ms2_dodd_lib', envir = environment())

                       cpd_obj <- cpd_dodd_lib
                       ms2_obj <- ms2_dodd_lib

                       rm(cpd_dodd_lib, ms2_dodd_lib);gc()
                     },
                     "stanford" = {
                       data('cpd_stanford_lib', envir = environment())
                       data("ms2_stanford_lib", envir = environment())
                       cpd_obj <- cpd_stanford_lib
                       ms2_obj <- ms2_stanford_lib
                       rm(cpd_stanford_lib, ms2_stanford_lib);gc()
                     },
                     "iroa" = {
                       data('cpd_iroa_lib', envir = environment())
                       data("ms2_iroa_lib", envir = environment())
                       cpd_obj <- cpd_iroa_lib
                       ms2_obj <- ms2_iroa_lib
                       rm(cpd_iroa_lib, ms2_iroa_lib);gc()
                     }
             )

             # load corresponding cpd information and ms2
             if (column == "hilic") {
               if (polarity == "positive") {
                 cpd_table <- cpd_obj$compound_table %>%
                   dplyr::select(lab_id:mz_neg, rt_hilic_pos, sources) %>%
                   dplyr::rename(rt = rt_hilic_pos)

                 if (lib == 'dodd') {
                   ms2_list <- ms2_obj$positive
                   db_info <- ms2_obj$lib_info
                 } else {
                   ms2_list <- ms2_obj$hilic_pos
                   db_info <- ms2_obj$lib_info
                 }

               } else {
                 cpd_table <- cpd_obj$compound_table %>%
                   dplyr::select(lab_id:mz_neg, rt_hilic_neg, sources) %>%
                   dplyr::rename(rt = rt_hilic_neg)

                 if (lib == 'dodd') {
                   ms2_list <- ms2_obj$negative
                   db_info <- ms2_obj$lib_info
                 } else {
                   ms2_list <- ms2_obj$hilic_neg
                   db_info <- ms2_obj$lib_info
                 }
               }
             } else {
               if (polarity == 'positive') {
                 cpd_table <- cpd_obj$compound_table %>%
                   dplyr::select(lab_id:mz_neg, rt_c18_pos, sources) %>%
                   dplyr::rename(rt = rt_c18_pos)

                 if (lib == 'dodd') {
                   ms2_list <- ms2_obj$positive
                   db_info <- ms2_obj$lib_info
                 } else {
                   ms2_list <- ms2_obj$c18_pos
                   db_info <- ms2_obj$lib_info
                 }
               } else {
                 cpd_table <- cpd_obj$compound_table %>%
                   dplyr::select(lab_id:mz_neg, rt_c18_neg, sources) %>%
                   dplyr::rename(rt = rt_c18_neg)

                 if (lib == 'dodd') {
                   ms2_list <- ms2_obj$negative
                   db_info <- ms2_obj$lib_info
                 } else {
                   ms2_list <- ms2_obj$c18_neg
                   db_info <- ms2_obj$lib_info
                 }
               }
             }

             # generate database info
             database.info <- list('Version' = db_info$version,
                                   "Source" = "MS",
                                   "Link" = "https://www.doddlab.org/",
                                   "Creater" = db_info$submitter,
                                   "Email" = "zhouzw@stanford.edu",
                                   "RT" = TRUE)

             # modify spectra.info
             spectra.info <- cpd_table %>%
               dplyr::mutate("Submitter" = db_info$submitter) %>%
               dplyr::rename("Lab.ID" = "lab_id",
                             "Compound.name" = "compound_name",
                             "mz" = "monoisotopic_mass",
                             "RT" = "rt",
                             "CAS.ID" = "cas_id",
                             "HMDB.ID" = "hmdb_id",
                             "KEGG.ID" = "kegg_id",
                             "Formula" = "formula",
                             "mz.pos" = "mz_pos",
                             "mz.neg" = "mz_neg") %>%
               dplyr::select("Lab.ID", "Compound.name", "mz", "RT", "CAS.ID",
                             "HMDB.ID", "KEGG.ID", "Formula", "mz.pos", "mz.neg",
                             "Submitter", dplyr::everything()) %>%
               dplyr::filter(!is.na(RT)) %>%
               dplyr::mutate(RT = round(RT*60, 2))


             # change ms2 list to data.frame
             ms2_list <- lapply(ms2_list, function(x){
               temp <- lapply(x, function(y){
                 as.data.frame(y)
               })

               names(temp) <- names(x)
               return(temp)
             })
             spectra.data <- list(Spectra.positive = ms2_list)


             db_output <- list('database.info' = database.info,
                               "spectra.info" = spectra.info,
                               "spectra.data" = spectra.data)

             message(crayon::green("Done!\n"))
             return(db_output)

           })


#' @title export_db_msdial
#' @author Zhiwei Zhou
#' @param lib dodd lab library. "stanford", "iroa"
#' @param column columne type. "hilic", "rp"
#' @param polarity ionization polarity. "positive", "negative"
#' @param msp_name ouput msp file name. e.g. "dodd_lib_stanford_hilic_pos.msp"
#' @param dir_path ouput directory path. Defaule: "."
#' @importFrom magrittr %>%
#' @importFrom crayon blue red yellow green bgRed
#' @importFrom stringr str_detect str_extract
#' @export
#' @examples
#' \dontrun{
#' export_db_msdial(lib = 'stanford',
#'                  column = 'hilic',
#'                  polarity = 'negative',
#'                  msp_name = 'stanford_hilic_neg.msp',
#'                  dir_path = '~/Project/00_IBD_project/Data/20230215_MS2_lib_processing/04_final_tables')
#' }

# export_db_msdial(lib = 'stanford',
#                  column = 'hilic',
#                  polarity = 'negative',
#                  msp_name = 'test.msp',
#                  dir_path = '~/Project/00_IBD_project/Data/20230215_MS2_lib_processing/04_final_tables')

setGeneric(name = 'export_db_msdial',
           def = function(lib = c('stanford', "iroa"),
                          column = c('hilic', 'c18'),
                          polarity = c('positive', 'negative'),
                          msp_name = "dodd_lib_stanford_hilic_pos.msp",
                          dir_path = '.'){

             message(crayon::blue("Start export database for MSDIAL...\n"))

             lib <- match.arg(lib)
             column <- match.arg(column)
             polarity <- match.arg(polarity)

             # load databases
             switch (lib,
                     # "doddLab" = {
                     #
                     # },
                     "stanford" = {
                       data('cpd_stanford_lib', envir = environment())
                       data("ms2_stanford_lib", envir = environment())
                       cpd_obj <- cpd_stanford_lib
                       ms2_obj <- ms2_stanford_lib
                       rm(cpd_stanford_lib, ms2_stanford_lib);gc()
                     },
                     "iroa" = {
                       data('cpd_iroa_lib', envir = environment())
                       data("ms2_iroa_lib", envir = environment())
                       cpd_obj <- cpd_iroa_lib
                       ms2_obj <- ms2_iroa_lib
                       rm(cpd_iroa_lib, ms2_iroa_lib);gc()
                     }
             )

             adduct_type <- ifelse(polarity == "positive", "[M+H]+", "[M-H]-")

             if (polarity == "positive" & column == "hilic") {
               temp_lib_id <- 'hilic_pos'
             } else if (polarity == "negative" & column == "hilic") {
               temp_lib_id <- 'hilic_neg'
             } else if (polarity == "positive" & column == "c18") {
               temp_lib_id <- 'c18_pos'
             } else {
               temp_lib_id <- 'c18_neg'
             }

             ms2_list <- ms2_obj[[temp_lib_id]]

             cpd_id_list <- names(ms2_list)
             walk(seq_along(cpd_id_list), function(i){
               temp_id <- cpd_id_list[i]
               idx <- which(cpd_obj$compound_table$lab_id == temp_id)
               temp_info <- cpd_obj$compound_table %>% dplyr::slice(i)

               temp_ms2_list <- ms2_list[[temp_id]]
               temp_ce_list <- names(temp_ms2_list)
               walk(seq_along(temp_ms2_list), function(i){
                 temp_spec <- temp_ms2_list[[i]]
                 temp_ce <- temp_ce_list[i]

                 generate_msp(file_name = file.path(dir_path, msp_name),
                              cmp_name = temp_info$compound_name,
                              precusormz = temp_info$mz_neg,
                              adduct = adduct_type,
                              instrument_type = "Q-TOF",
                              instrument = "Agilent-6545-XF",
                              smiles = temp_info$smiles,
                              inchikey = temp_info$inchikey,
                              formula = temp_info$formula,
                              polarity = polarity,
                              ce = temp_ce,
                              rt = ifelse(is.na(temp_info$rt_c18_pos), "", temp_info$rt_c18_pos),
                              doddlib_id = temp_info$lab_id,
                              kegg_id = temp_info$kegg_id,
                              hmdb_id = temp_info$hmdb_id,
                              pubchem_id = temp_info$pubchem_id,
                              spec = temp_spec)

               })
             })

             message(crayon::green("Done!\n"))

           })






#' @title generate_msp
#' @author Zhiwei Zhou
#' @description Convert MS/MS spectra to MSP files
#' @param file_name Required. The file name of MSP. The suffix of file name must be ".msp".
#' @param cmp_name Required
#' @param precusormz Required
#' @param spec Required. It should be in a dataframe, 1st: "mz", 2nd: "intensity"
#' @param adduct Default: NULL
#' @param instrument_type Default: NULL
#' @param instrument Default: NULL
#' @param smiles Default: NULL
#' @param inchikey Default: NULL
#' @param inchikey1 Default: NULL
#' @param formula Default: NULL
#' @param polarity 'Positive' or 'Negative'
#' @param ce Default: NULL
#' @param rt Default: NULL
#' @param ccs Default: NULL
#' @param doddlib_id Default: NULL
#' @param kegg_id Default: NULL
#' @param hmdb_id Default: NULL
#' @param pubchem_id Default: NULL
#' @param links Default: ''
#' @param comment Default: ''
#' @export
#' @examples
#' setwd('F:/01 MetIMMS/00 data processing/190515 external validation msms data extraction/')
#'     GenerateMSP(file_name = 'zhumetlib_validation_pos_20v_190520.msp',
#'                 cmp_name = external_data_pos$compound_name[i],
#'                 precusormz = external_data_pos$mz[i],
#'                 adduct = external_data_pos$adducts[i],
#'                 instrument_type = 'LC-ESI-qTof',
#'                 instrument = 'LC-ESI-qTof',
#'                 smiles = external_data_pos$smiles[i],
#'                 inchikey = external_data_pos$inchikey[i],
#'                 inchikey1 = external_data_pos$inchikey1[i],
#'                 formula = external_data_pos$formula[i],
#'                 polarity = 'Positive',
#'                 ce = '20',
#'                 ccs = external_data_pos$CCS[i],
#'                 zhulib_id = external_data_pos$matched_zhulib_id[i],
#'                 pubchem_id = external_data_pos$pubchem_cid[i],
#'                 comment = paste(external_data_pos$id[i], external_data_pos$source[i], sep = ' '),
#'                 spec = temp_spec)



setGeneric(name = 'generate_msp',
           def = function(
    file_name = "./test.msp",
    cmp_name = NULL,
    precusormz = NULL,
    adduct=NULL,
    instrument_type=NULL,
    instrument=NULL,
    smiles=NULL,
    inchikey=NULL,
    inchikey1=NULL,
    formula=NULL,
    polarity=c('positive', 'negative'),
    ce=NULL,
    rt=NULL,
    ccs=NULL,
    doddlib_id=NULL,
    kegg_id=NULL,
    hmdb_id=NULL,
    pubchem_id=NULL,
    links='',
    comment='',
    spec=NULL
           ){

             if (!stringr::str_detect(file_name, '.msp')) {
               stop('The suffix of file name must be .msp\n')
             }

             if (is.null(cmp_name)) {
               stop('Please input cmp_name\n')
             }

             if (is.null(precusormz)) {
               stop('Please input precusormz\n')
             }

             if (is.null(spec)) {
               stop('Please input spec\n')
             }

             polarity = match.arg(polarity)


             if (ncol(spec)!=2) {
               stop('Please check the format of spectrum. It should be in a dataframe, 1st: "mz", 2nd: "intensity" \n')
             }

             if (!all(colnames(spec)==c('mz', 'intensity'))) {
               stop('Please check the format of spectrum. It should be in a dataframe, 1st: "mz", 2nd: "intensity" \n')
             }



             # write into msp
             file_result <- file(description = file_name, open = "a")
             cat('NAME: ', cmp_name, '\n', sep = '', file = file_result)
             cat('PRECURSORMZ: ', precusormz, '\n', sep = '', file = file_result)

             if (!is.null(adduct)) {
               cat('PRECURSORTYPE: ', adduct, '\n', sep = '', file = file_result)
             }

             if (!is.null(instrument_type)) {
               cat('INSTRUMENTTYPE: ', instrument_type, '\n', sep = '', file = file_result)
             }

             if (!is.null(instrument)) {
               cat('INSTRUMENT: ', instrument, '\n', sep = '', file = file_result)
             }

             if (!is.null(smiles)) {
               cat('SMILES: ', smiles, '\n', sep = '', file = file_result)
             }

             if (!is.null(inchikey)) {
               cat('INCHIKEY: ', inchikey, '\n', sep = '', file = file_result)
             }

             if (!is.null(inchikey1)) {
               cat('INCHIKEY1: ', inchikey1, '\n', sep = '', file = file_result)
             }

             if (!is.null(formula)) {
               cat('FORMULA: ', formula, '\n', sep = '', file = file_result)
             }

             if (!is.null(polarity)) {
               cat('IONMODE: ', polarity, '\n', sep = '', file = file_result)
             }

             if (!is.null(ce)) {
               cat('COLLISIONENERGY: ', ce, '\n', sep = '', file = file_result)
             }

             if (!is.null(rt)) {
               cat('RETENTIONTIME: ', rt, '\n', sep = '', file = file_result)
             }

             if (!is.null(ccs)) {
               cat('COLLISIONCROSSSECTION: ', ccs, '\n', sep = '', file = file_result)
             }

             if (!is.null(doddlib_id)) {
               cat('DODDLAB: ', doddlib_id, '\n', sep = '', file = file_result)
             }

             if (!is.null(kegg_id)) {
               cat('KEGG: ', kegg_id, '\n', sep = '', file = file_result)
             }

             if (!is.null(hmdb_id)) {
               cat('HMDB: ', hmdb_id, '\n', sep = '', file = file_result)
             }

             if (!is.null(pubchem_id)) {
               cat('PUBCHEM: ', pubchem_id, '\n', sep = '', file = file_result)
             }

             cat('Links: ', links, '\n', sep = '', file = file_result)
             cat('Comment: ', comment, "\n", sep = '', file = file_result)

             cat('Num Peaks: ',  nrow(spec),  '\n',  sep = '', file = file_result)

             for (i in 1:nrow(spec)) {
               cat(paste(as.numeric(round(spec[i,1], digits = 4)),
                         as.numeric(round(spec[i,2], digits = 2)),
                         collapse = ' '),
                   '\n', sep = '', file = file_result)
             }

             cat('\n', file = file_result)

             close(file_result)
           }
)




#' @title export_ms2_db
#' @author Zhiwei Zhou
#' @param lib dodd lab library. "dodd"
#' @param polarity ionization polarity. "positive", "negative"
#' @param ms1_version MS1 version. "v2.6.1", "v2.5.0", "v2.4.0", "v2.0.0"
#' @param ms2_version MS2 version. "v2.6.1", "v2.5.0", "v2.4.0", "v2.3.0"
#' @param msp_name ouput msp file name. e.g. "dodd_lib_stanford_hilic_pos.msp"
#' @param dir_path ouput directory path. Defaule: "."
#' @importFrom magrittr %>%
#' @importFrom crayon blue red yellow green bgRed
#' @importFrom stringr str_detect str_extract
#' @export
#' @examples
#' \dontrun{
#' export_ms2_db(lib = 'dodd',
#'               ms1_version = 'v2.6.1',
#'               ms2_version = 'v2.6.1',
#'               polarity = 'positive',
#'               msp_name = "dodd_lib_pos_v2.6.1.msp",
#'               dir_path = '~/Project/04_package/00_Database/DoddLib/06_exported_db/')
#'
#' export_ms2_db(lib = 'dodd',
#'               ms1_version = 'v2.6.1',
#'               ms2_version = 'v2.6.1',
#'               polarity = 'negative',
#'               msp_name = "dodd_lib_neg_v2.6.1.msp",
#'               dir_path = '~/Project/04_package/00_Database/DoddLib/06_exported_db/')
#' }


# lib <- 'dodd'
# ms1_version <- 'v2.6.1'
# ms2_version <- 'v2.6.1'
# polarity <- 'positive'
# msp_name <- "dodd_lib_pos.msp"
# dir_path <- '~/Project/04_package/00_Database/DoddLib/06_exported_db/'

# export_ms2_db(lib = 'dodd',
#               ms1_version = 'v2.6.1',
#               ms2_version = 'v2.6.1',
#               polarity = 'positive',
#               msp_name = "dodd_lib_pos_v2.6.1.msp",
#               dir_path = '~/Project/04_package/00_Database/DoddLib/06_exported_db/')
#
# export_ms2_db(lib = 'dodd',
#               ms1_version = 'v2.6.1',
#               ms2_version = 'v2.6.1',
#               polarity = 'negative',
#               msp_name = "dodd_lib_neg_v2.6.1.msp",
#               dir_path = '~/Project/04_package/00_Database/DoddLib/06_exported_db/')

setGeneric(name = 'export_ms2_db',
           def = function(lib = c('dodd'),
                          ms1_version = c('v2.6.1', 'v2.5.0', 'v2.4.0', 'v2.0.0'),
                          ms2_version = c('v2.6.1', 'v2.5.0', 'v2.4.0', 'v2.3.0'),
                          polarity = c('positive', 'negative'),
                          msp_name = "dodd_lib_pos.msp",
                          dir_path = '.'){

             message(crayon::blue("Start export database for MSDIAL...\n"))

             lib <- match.arg(lib)
             polarity <- match.arg(polarity)
             ms1_version <- match.arg(ms1_version)
             ms2_version <- match.arg(ms2_version)

             # load databases
             data('list_cpd_dodd_lib', envir = environment())
             cpd_obj <- list_cpd_dodd_lib[[ms1_version]]
             rm(list_cpd_dodd_lib);gc()

             data('list_ms2_dodd_lib', envir = environment())
             ms2_obj <- list_ms2_dodd_lib[[ms2_version]][[polarity]]
             rm(list_ms2_dodd_lib);gc()

             adduct_type <- ifelse(polarity == "positive", "[M+H]+", "[M-H]-")

             lab_id_ms2 <- names(ms2_obj)
             walk(seq_along(lab_id_ms2), function(i){
               temp_id <- lab_id_ms2[i]
               idx <- which(cpd_obj$compound_table$lab_id == temp_id)
               temp_info <- cpd_obj$compound_table %>% dplyr::slice(idx)

               temp_ms2_list <- ms2_obj[[temp_id]]
               temp_ce_list <- names(temp_ms2_list)
               walk(seq_along(temp_ms2_list), function(i){
                 temp_spec <- temp_ms2_list[[i]]
                 temp_ce <- temp_ce_list[i]

                 generate_msp(file_name = file.path(dir_path, msp_name),
                              cmp_name = temp_info$compound_name,
                              precusormz = ifelse(polarity == "positive", temp_info$mz_pos, temp_info$mz_neg),
                              adduct = adduct_type,
                              instrument_type = "Q-TOF",
                              instrument = "Agilent-6545-XF",
                              smiles = temp_info$smiles,
                              inchikey = temp_info$inchikey,
                              formula = temp_info$formula,
                              polarity = polarity,
                              ce = temp_ce,
                              # rt = ifelse(is.na(temp_info$rt_c18_pos), "", temp_info$rt_c18_pos),
                              doddlib_id = temp_info$lab_id,
                              kegg_id = temp_info$kegg_id,
                              hmdb_id = temp_info$hmdb_id,
                              pubchem_id = temp_info$pubchem_id,
                              spec = temp_spec)

               })
             })

             message(crayon::green("Done!\n"))
           })
