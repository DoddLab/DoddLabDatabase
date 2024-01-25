################################################################################
# retrieve_raw_ms2 -------------------------------------------------------------
#' @title retrieve_raw_ms2
#' @author Zhiwei Zhou
#' @param file mzml or mzxml file. ms2 files
#' @param isolation_window 'narrow', 'medium', 'wide'
#' @param dir_path '.'
#' @param is_purification whether do the ms2 purification
#' @importFrom magrittr %>%
#' @importFrom crayon blue red yellow green bgRed
#' @importFrom stringr str_detect str_extract
#' @export
#' @examples
#' \dontrun{
#' file <- '~/Project/00_IBD_project/Data/20230428_Zhiwei_StanfordPool_MS2/mzxml_haoqing_data/Haoqing_MS2_data/MS2_lib_mzXML/Stanford_HILIC_pos/SU-A01_1.mzXML'
#' file <- '~/Project/00_IBD_project/Data/20230428_Zhiwei_StanfordPool_MS2/mzml_hilic_pos/StanfordPool1-2.mzML'
#'
#' test <- retrieve_raw_ms2(file = 'SU-A01_1.mzXML',
#'                          isolation_window = 'narrow',
#'                          dir_path = '~/Project/00_IBD_project/Data/20230428_Zhiwei_StanfordPool_MS2/mzxml_haoqing_data/Haoqing_MS2_data/MS2_lib_mzXML/Stanford_HILIC_pos')
#' }

# file <- '~/Project/00_IBD_project/Data/20230428_Zhiwei_StanfordPool_MS2/mzxml_haoqing_data/Haoqing_MS2_data/MS2_lib_mzXML/Stanford_HILIC_pos/SU-A01_1.mzXML'
# file <- '~/Project/00_IBD_project/Data/20230428_Zhiwei_StanfordPool_MS2/mzml_hilic_pos/StanfordPool1-2.mzML'
#
# test <- retrieve_raw_ms2(file = 'SU-A01_1.mzXML',
#                          isolation_window = 'narrow',
#                          dir_path = '~/Project/00_IBD_project/Data/20230428_Zhiwei_StanfordPool_MS2/mzxml_haoqing_data/Haoqing_MS2_data/MS2_lib_mzXML/Stanford_HILIC_pos')
# test2 <- retrieve_raw_ms2(file = 'SU-A01_1.mzXML',
#                           isolation_window = 'narrow',
#                           dir_path = '~/Project/00_IBD_project/Data/20230428_Zhiwei_StanfordPool_MS2/mzxml_haoqing_data/Haoqing_MS2_data/MS2_lib_mzXML/Stanford_HILIC_pos',
#                           is_purification = FALSE)
#
#
# test <- retrieve_raw_ms2(file = 'StanfordPool1-2.mzML',
#                          isolation_window = 'narrow',
#                          dir_path = '~/Project/00_IBD_project/Data/20230428_Zhiwei_StanfordPool_MS2/mzml_hilic_pos')
#
# test2 <- retrieve_raw_ms2(file = 'StanfordPool1-2.mzML',
#                          isolation_window = 'narrow',
#                          dir_path = '~/Project/00_IBD_project/Data/20230428_Zhiwei_StanfordPool_MS2/mzml_hilic_pos',
#                          is_purification = FALSE)

retrieve_raw_ms2 <- function(file,
                             isolation_window = c('narrow', 'medium', 'wide'),
                             dir_path = '.',
                             is_purification = TRUE) {
  isolation_window <- match.arg(isolation_window)
  switch (isolation_window,
    'narrow' = {
      isolation_range <- 0.65
    },
    'medium' = {
      isolation_range <- 4
    },
    'wide' = {
      isolation_range <- 9
    }
  )

  mzml.data <- mzR::openMSfile(file.path(dir_path, file))
  mzml.info <- mzR::header(mzml.data)
  mzml.peak <- mzR::peaks(mzml.data)

  mzml_info_ms1 <- mzml.info %>%
    dplyr::filter(msLevel == 1)
  mzml_info_ms2 <- mzml.info %>%
    dplyr::filter(msLevel == 2)


  # In targeted ms/ms, no precursorScanNum provided. Assign the closed ms1 scan as precursor scan
  precursor_scan_num <- mzml_info_ms2$precursorScanNum
  if (all(is.na(precursor_scan_num))) {
    temp_position <- outer(mzml_info_ms1$seqNum, mzml_info_ms2$seqNum, function(x, y){
      x - y
    })
    # the ms1 scan number must less than ms2 scan. assign 10000 if the ms1 scan larger than ms2 scan
    precursor_idx <- apply(temp_position, 2, function(x){
      x[x>0] <- 10000
      which.min(abs(x))
    })

    precursor_seq_num <- mzml_info_ms1$seqNum[precursor_idx]
    precursor_rt <- mzml_info_ms1$retentionTime[precursor_idx]
  } else {
    precursor_seq_num <- match(precursor_scan_num, mzml_info_ms1$acquisitionNum) %>% mzml_info_ms1$seqNum[.]
    precursor_rt <- match(precursor_scan_num, mzml_info_ms1$acquisitionNum) %>% mzml_info_ms1$retentionTime[.]
   }


  precursor_info <- mapply(function(seq_num, mz_precursor){
    temp_spec <- mzml.peak[[seq_num]]

    temp <- which(abs(temp_spec[,'mz'] - mz_precursor) <= 0.0015)
    if (length(temp) < 1) {
      precursor_info <- tibble::tibble(precursor_mz = mz_precursor,
                                       intensity = 0,
                                       purity = -1,
                                       ms1_scan = seq_num,
                                       mz_error = -1)

      return(precursor_info)
    }

    # select the correct precursor
    idx <- which.min(abs(temp_spec[,'mz'] - mz_precursor))
    int_precursor <- temp_spec[,'intensity'][idx]
    mz_exp <- temp_spec[,'mz'][idx]
    mz_error <- abs(temp_spec[,'mz'] - mz_precursor)[idx]

    # calculate the precursor purity
    temp_idx <- which((temp_spec[,'mz'] >= (mz_precursor - isolation_range)) &
                        (temp_spec[,'mz'] <= (mz_precursor + isolation_range)))
    purity_precursor <- int_precursor/sum(temp_spec[,'intensity'][temp_idx])

    # generate precursor info
    precursor_info <- tibble::tibble(precursor_mz = mz_precursor,
                                     intensity = int_precursor,
                                     purity = purity_precursor,
                                     ms1_scan = seq_num,
                                     mz_error = mz_error)

    return(precursor_info)
  },
  seq_num = precursor_seq_num,
  mz_precursor = mzml_info_ms2$precursorMZ,
  SIMPLIFY = FALSE)

  precursor_info <- precursor_info %>% dplyr::bind_rows()

  ms2_info <- mzml_info_ms2 %>%
    dplyr::bind_cols(precursor_info) %>%
    dplyr::mutate(precursor_rt = precursor_rt) %>%
    dplyr::mutate(precursor_mz = precursor_mz,
                  ms2_rt = retentionTime,
                  ce = collisionEnergy,
                  ms2_scan = seqNum,
                  precuror_int = intensity,
                  ms2_file = file) %>%
    dplyr::select(precursor_mz,
                  precursor_rt,
                  precuror_int,
                  ms2_rt,
                  ce,
                  purity,
                  ms1_scan,
                  ms2_scan,
                  ms2_file)

  ms2_spec <- mzml.peak[ms2_info$ms2_scan]

  if (is_purification) {
    purified_spec <- mapply(function(x, y){
      DoddLabMetID::purifyMs2(x,
                              mz_precursor = y,
                              is_include_precursor = TRUE,
                              is_remove_ring_effect = TRUE,
                              ppm_precursor_filter = 10,
                              int_ms2_min_abs = 50,
                              int_ms2_min_relative = 0.01)
    },
    x = ms2_spec,
    y = ms2_info$precursor_mz,
    SIMPLIFY = FALSE)

    idx_eff <- which(sapply(purified_spec, length) > 0)
    if (length(idx_eff) > 0) {
      purified_spec <- purified_spec[idx_eff]
      ms2_info <- ms2_info[idx_eff,]
      result <- list('info' = ms2_info,
                     'spec' = purified_spec)
      return(result)
    } else{
      result <- list(spec_info = NULL,
                     best_spec = NULL)
      return(result)
    }
  }

  result <- list('info' = ms2_info,
                 'spec' = ms2_spec)

}



# group_ms2 --------------------------------------------------------------------
#' @title group_ms2
#' @author Zhiwei Zhou
#' @param target_list a target list (data.frame). First 3 columns, id, mz, rt
#' @param ms2_data ms2 data. Generated by retrieve_ms2
#' @param label Default: ''. Usually can label batch name, e.g. '2023_Zhiwei'
#' @param mode Default: ''. Usually can label column and polarity, e.g. 'hilic_pos'
#' @param mz_tol Default: 15. Match target list with ms2_info tolerance
#' @param rt_tol Default: 15. Match target list with ms2_info tolerance
#' @importFrom magrittr %>%
#' @importFrom crayon blue red yellow green bgRed
#' @importFrom stringr str_detect str_extract
#' @export
#' @examples
#' \dontrun{
#' }

# ms2_data <- retrieve_raw_ms2(file = 'SU-A01_1.mzXML',
#                              isolation_window = 'narrow',
#                              dir_path = '~/Project/00_IBD_project/Data/20230428_Zhiwei_StanfordPool_MS2/mzxml_haoqing_data/Haoqing_MS2_data/MS2_lib_mzXML/Stanford_HILIC_pos')
#
# generate target list
# load('~/Project/00_IBD_project/Data/20230428_Zhiwei_StanfordPool_MS2/mzxml_haoqing_data/Haoqing_MS2_data/MS2_lib_processing/01_compound_info/ddodd_stanford_lib_hilic_pos_230209.RData')
# mix_id <- 'SU-A01_1.mzXML' %>% stringr::str_extract('[A-Z]\\d+')
# target_list <- ddodd_stanford_lib_hilic_pos %>%
#   dplyr::filter(Stanford_Pool == mix_id) %>%
#   dplyr::filter(!is.na(RT)) %>%
#   dplyr::select(Lab.ID, mz.pos, RT) %>%
#   dplyr::rename(id = Lab.ID,
#                 mz = mz.pos,
#                 rt = RT)
# save(target_list,
#      file = '~/Project/00_IBD_project/Data/20230428_Zhiwei_StanfordPool_MS2/mzxml_haoqing_data/target_list')
#
# test <- group_ms2(target_list = target_list,
#                   ms2_data = ms2_data,
#                   label = '2021_Haoqing',
#                   mode = 'hilic_pos',
#                   mz_tol = 15,
#                   rt_tol = 15)
#
# ms2_data <- retrieve_raw_ms2(file = 'StanfordPool1-2.mzML',
#                              isolation_window = 'narrow',
#                              dir_path = '~/Project/00_IBD_project/Data/20230428_Zhiwei_StanfordPool_MS2/mzml_hilic_pos')
#
# low_quality_ms2_list <- readxl::read_xlsx('~/Project/00_IBD_project/Data/20230428_Zhiwei_StanfordPool_MS2/low_quality_ms2_std_stanford_230412.xlsx')
# low_quality_ms2_list <- low_quality_ms2_list %>%
#   dplyr::mutate(pool_sample_id = stringr::str_replace(pool_sample_id, pattern = 'StandfordPool', replacement = 'StanfordPool'))
# target_list <- low_quality_ms2_list %>%
#   dplyr::filter(pool_sample_id == 'StanfordPool1') %>%
#   dplyr::filter(!is.na(rt_hilic_pos)) %>%
#   dplyr::mutate(rt = rt_hilic_pos*60) %>%
#   dplyr::select(lab_id, mz_pos, rt) %>%
#   dplyr::rename(id = lab_id,
#                 mz = mz_pos,
#                 rt = rt)
#
# test <- group_ms2(target_list = target_list,
#                   ms2_data = ms2_data,
#                   label = '2023_zhiwei',
#                   mode = 'hilic_pos',
#                   mz_tol = 15,
#                   rt_tol = 15)

group_ms2 <- function(target_list,
                      ms2_data,
                      label = '',
                      mode = 'hilic_pos',
                      mz_tol = 15,
                      rt_tol = 15,
                      rt_error_type = 'abs'){
  # browser()
  if (missing(ms2_data)) {
    stop('Plese input the ms2_data. It should be run first')
  }

  warning <- check_conflict(target_list = target_list)
  target_list <- target_list %>%
    dplyr::mutate(warning = warning)

  temp_info <- target_list %>% dplyr::select(mz, rt)
  temp_ms2_info <- ms2_data$info %>% dplyr::select(precursor_mz, precursor_rt)
  # match acquired ms2 with target list
  idx_table <- masstools::mz_rt_match(data1 = as.data.frame(temp_info),
                                      data2 = as.data.frame(temp_ms2_info),
                                      mz.tol = mz_tol,
                                      rt.tol = rt_tol,
                                      rt.error.type = rt_error_type)

  if (length(idx_table) == 0) {
    return(NULL)
  }

  # split spectra into list, group by each standard
  temp_idx <- match(idx_table$Index1, unique(idx_table$Index1))
  idx_corresponding <- split(idx_table, temp_idx) %>% as.list()
  id_list <- unique(idx_table$Index1) %>% target_list$id[.]

  target_mz_list <- unique(idx_table$Index1) %>% target_list$mz[.]
  target_rt_list <- unique(idx_table$Index1) %>% target_list$rt[.]
  warning_list <- unique(idx_table$Index1) %>% target_list$warning[.]

  result <- lapply(seq_along(idx_corresponding), function(i){
    temp_idx <- idx_corresponding[[i]]$Index2
    new_ms2_info <- ms2_data$info[temp_idx,] %>%
      dplyr::mutate(id = id_list[i],
                    target_mz = target_mz_list[i],
                    target_rt = target_rt_list[i],
                    label = label,
                    mode = mode,
                    warning = warning_list[i]) %>%
      dplyr::mutate(mz_error = abs(target_mz - precursor_mz),
                    rt_error = abs(target_rt - precursor_rt)) %>%
      dplyr::select(id, target_mz, target_rt, precursor_mz:purity, mz_error, rt_error,
                    ms1_scan:ms2_scan, label, mode, warning, ms2_file)

    new_spec <- ms2_data$spec[temp_idx]

    result <- list(info = new_ms2_info,
                   spec = new_spec)

    return(result)
  })

  names(result) <- id_list

  return(result)
}



check_conflict <- function(target_list,
                           mz_tol_conflict = 0.65,
                           rt_tol_conflict = 15) {
  # browser()
  mz_error_matrix <- outer(target_list$mz, target_list$mz, function(x, y){
    abs(x - y)
  })

  rt_error_matrix <- outer(target_list$rt, target_list$rt, function(x, y){
    abs(x - y)
  })

  idx_mz_error_matrix <- which(mz_error_matrix <= mz_tol_conflict, arr.ind = TRUE)
  idx_rt_error_matrix <- which(rt_error_matrix <= rt_tol_conflict, arr.ind = TRUE)

  idx_mz_conflict <- which(!(idx_mz_error_matrix[,'row'] == idx_mz_error_matrix[,'col']))
  idx_rt_conflict <- which(!(idx_rt_error_matrix[,'row'] == idx_rt_error_matrix[,'col']))

  label <- rep('', length(target_list$mz))
  if (length(idx_mz_conflict) > 0 & length(idx_rt_conflict) > 0) {
    temp <- intersect(idx_mz_error_matrix[idx_mz_conflict], idx_rt_error_matrix[idx_rt_conflict])
    if (length(temp) > 0) {
      label[temp] <- 'conflict'
    }
  }

  return(label)
}


################################################################################

#' @title merge_ms2
#' @author Zhiwei Zhou
#' @description merge the ms2 list, and organized with ID
#' @param ms2_list a list of STD ms2 data.
#' @importFrom magrittr %>%
#' @importFrom crayon blue red yellow green bgRed
#' @importFrom stringr str_detect str_extract
#' @export
#' @examples
#' \dontrun{
#' load('~/Project/00_IBD_project/Data/20230428_Zhiwei_StanfordPool_MS2/mzxml_haoqing_data/Haoqing_MS2_data/MS2_lib_230509/ms2_stanford_C18_pos_raw_230509.RData')
#' load('~/Project/00_IBD_project/Data/20230428_Zhiwei_StanfordPool_MS2/mzxml_haoqing_data/Haoqing_MS2_data/MS2_lib_230509/ms2_stanford_HILIC_pos_raw_230509.RData')
#'
#' # combine these two ms2 list as one list
#' ms2_list <- c(ms2_stanford_HILIC_pos_raw, ms2_stanford_C18_pos_raw)
#' ms2_list <- merge_ms2(ms2_list)
#' }


# load('~/Project/00_IBD_project/Data/20230428_Zhiwei_StanfordPool_MS2/mzxml_haoqing_data/Haoqing_MS2_data/MS2_lib_230509/ms2_stanford_C18_pos_raw_230509.RData')
# load('~/Project/00_IBD_project/Data/20230428_Zhiwei_StanfordPool_MS2/mzxml_haoqing_data/Haoqing_MS2_data/MS2_lib_230509/ms2_stanford_HILIC_pos_raw_230509.RData')
# ms2_list <- c(ms2_stanford_HILIC_pos_raw, ms2_stanford_C18_pos_raw)
# ms2_list <- merge_ms2(ms2_list)

merge_ms2 <- function(ms2_list) {
  id_list <- lapply(ms2_list, function(x){
    names(x)
  })

  unique_id <- id_list %>% unlist() %>% unique() %>% sort()

  result <- pbapply::pblapply(unique_id, function(temp_id){
    idx <- lapply(id_list, function(x){
      which(x == temp_id)
    })

    temp_data <- mapply(function(x, y){
      if (length(x) == 0) return(NULL)
      return(x[y])
    },
    x = ms2_list,
    y = idx,
    SIMPLIFY = FALSE)

    info1 <- lapply(temp_data, function(x){
      if (length(x) == 0) {
        return(NULL)
      }

      x[[1]]$info
    }) %>%
      dplyr::bind_rows()

    spec1 <- lapply(temp_data, function(x){
      if (length(x) == 0) {
        return(NULL)
      }

      x[[1]]$spec
    }) %>% do.call(c, .)

    result <- list(info = info1,
                   spec = spec1)

    return(result)
  })

  names(result) <- unique_id
  return(result)
}



#' @title select_best_ms2
#' @author Zhiwei Zhou
#' @param ms2_list a merged ms2 list. The merged ms2 list can be generated by merge_ms2()
#' @param precursor_mz_error_cutoff the precursor mz error cutoff. Default: 0.0015 Da
#' @param precursor_rt_error_cutoff the precursor rt error cutoff. Default: 15 second
#' @param precursor_int_cutoff the minimum precursor intensity. Default: 500 counts
#' @param precursor_purity_cutoff the minmium precursor purity. Default: 0.9
#' @param is_keep_conflict_std whether keep the possible conflict ms2 spec. e.g. some targeted isobaric/isomers run in one injection, and have very similar RT. Default: TRUE
#' @param rank_by the way to select the best spec, includes "intensity", "intensity_purity", "rt_error". Default: "intensity"
#' @importFrom magrittr %>%
#' @importFrom crayon blue red yellow green bgRed
#' @importFrom stringr str_detect str_extract
#' @export
#' @examples
#' \dontrun{
#' load('~/Project/00_IBD_project/Data/20230428_Zhiwei_StanfordPool_MS2/mzxml_haoqing_data/Haoqing_MS2_data/MS2_lib_230509/ms2_stanford_C18_pos_raw_230509.RData')
#' load('~/Project/00_IBD_project/Data/20230428_Zhiwei_StanfordPool_MS2/mzxml_haoqing_data/Haoqing_MS2_data/MS2_lib_230509/ms2_stanford_HILIC_pos_raw_230509.RData')
#'
#' # combine these two ms2 list as one list
#' ms2_list <- c(ms2_stanford_HILIC_pos_raw, ms2_stanford_C18_pos_raw)
#' ms2_list <- merge_ms2(ms2_list)
#' select_best_ms2(ms2_list = ms2_list,
#'                 precursor_mz_error_cutoff = 0.0015, # absolute mz error, Da
#'                 precursor_rt_error_cutoff = 15, # absolute rt error, second
#'                 precursor_int_cutoff = 500,
#'                 precursor_purity_cutoff = 0.9,
#'                 is_keep_conflict_std = TRUE,
#'                 rank_by = "intensity_purity")
#' }

# temp_ms2_list <- ms2_list[[1]]
# select_best_ms2(ms2_list = ms2_list,
#                 precursor_mz_error_cutoff = 0.0015, # absolute mz error, Da
#                 precursor_rt_error_cutoff = 15, # absolute rt error, second
#                 precursor_int_cutoff = 500,
#                 precursor_purity_cutoff = 0.9,
#                 is_keep_conflict_std = TRUE,
#                 rank_by = "intensity_purity")

select_best_ms2 <- function(ms2_list,
                            precursor_mz_error_cutoff = 0.0015, # absolute mz error, Da
                            precursor_rt_error_cutoff = 15, # absolute rt error, second
                            precursor_int_cutoff = 500,
                            precursor_purity_cutoff = 0.9,
                            is_keep_conflict_std = TRUE,
                            rank_by = c("intensity", "intensity_purity", "rt_error")){

  rank_by <- match.arg(rank_by)

  best_ms2_list <- pbapply::pblapply(ms2_list, function(temp_ms2_list){
    temp_ms2_info <- temp_ms2_list$info %>%
      dplyr::mutate(index = seq(dplyr::n())) %>%
      dplyr::mutate(intensity_purity = precuror_int*purity)

    temp_ms2_info <- temp_ms2_info %>%
      dplyr::filter(mz_error <= precursor_mz_error_cutoff) %>%
      dplyr::filter(rt_error <= precursor_rt_error_cutoff) %>%
      dplyr::filter(precuror_int >= precursor_int_cutoff) %>%
      dplyr::filter(purity >= precursor_purity_cutoff)

    if (!is_keep_conflict_std) {
      temp_ms2_info <- temp_ms2_info %>%
        dplyr::filter(warning != 'conflict')
    }

    if (nrow(temp_ms2_info) < 1) {
      return(NULL)
    }

    ces <- unique(temp_ms2_info$ce) %>% sort()
    # for multiple spec for same ce, select the spec with the highest precursor intensity:
    best_ms2_info <- lapply(ces, function(temp_ce){
      idx <- which(temp_ms2_info$ce == temp_ce)
      temp_info <- temp_ms2_info[idx,]

      switch (rank_by,
        "intensity" = {
          temp_info <- temp_info %>%
            dplyr::arrange(dplyr::desc(descprecuror_int)) %>%
            dplyr::slice(1)
        },
        "intensity_purity" = {
          temp_info <- temp_info %>%
            dplyr::arrange(dplyr::desc(intensity_purity)) %>%
            dplyr::slice(1)
        },
        "rt_error" = {
          temp_info <- temp_info %>%
            dplyr::arrange(rt_error) %>%
            dplyr::slice(1)
        }
      )

      return(temp_info)
    })

    best_ms2_info <- best_ms2_info %>% dplyr::bind_rows()
    temp_ms2_idx <- best_ms2_info$index
    best_ms2_spec <- temp_ms2_list$spec[temp_ms2_idx]
    names(best_ms2_spec) <- best_ms2_info$ce

    result <- list(info = best_ms2_info,
                   spec = best_ms2_spec)

    return(result)
  })

  return(best_ms2_list)
}
