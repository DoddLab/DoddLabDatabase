# ################################################################################
# # compound db dodd lib ---------------------------------------------------------
# load('~/Project/04_package/00_Database/DoddLib/01_compound_RT_lib/02_version_230221/cpd_dodd_lib_220221.RData')
# usethis::use_data(cpd_dodd_lib, overwrite = TRUE)
#
# load('~/Project/04_package/00_Database/DoddLib/01_compound_RT_lib/02_version_230221/ms2_dodd_lib_220221.RData')
# usethis::use_data(ms2_dodd_lib, overwrite = TRUE)
#
# ################################################################################
# # separated dbs ----------------------------------------------------------------
# load('~/Project/04_package/00_Database/DoddLib/01_compound_RT_lib/02_version_230221/compound_lib_stanford_230221.RData')
# usethis::use_data(cpd_stanford_lib, overwrite = TRUE)
#
# load('~/Project/04_package/00_Database/DoddLib/01_compound_RT_lib/02_version_230221/compound_lib_iroa_230221.RData')
# usethis::use_data(cpd_iroa_lib, overwrite = TRUE)
#
# load('~/Project/04_package/00_Database/DoddLib/02_ms2_lib/v230216/ms2_stanford_lib_230215.RData')
# usethis::use_data(ms2_stanford_lib, overwrite = TRUE)
#
# load('~/Project/04_package/00_Database/DoddLib/02_ms2_lib/v230216/ms2_iroa_lib_230215.RData')
# usethis::use_data(ms2_iroa_lib, overwrite = TRUE)
#
# ################################################################################
# # internal standards -----------------------------------------------------------
# load('~/Project/04_package/00_Database/DoddLib/03_isotope_standard_ISTD/v230215/istd_lib_230215.RData')
# usethis::use_data(istd_lib, overwrite = TRUE)
#
# ################################################################################
# MS-DIAL ----------------------------------------------------------------------
# load('~/Project/04_package/00_Database/MSDIAL_library/v230302/msdial_lib_230306.RData')
# usethis::use_data(msdial_lib, overwrite = TRUE)E
#
# ################################################################################
# adduct_nl lib ----------------------------------------------------------------
# load('~/Project/04_package/MetDNA2/data/lib_adduct_nl.rda')
# usethis::use_data(lib_adduct_nl, overwrite = TRUE)
#
# ################################################################################
# ms2_dodd_lib -----------------------------------------------------------------
# load('~/Project/04_package/00_Database/DoddLib/02_ms2_lib/v230511_v2.1.0/ms2_dodd_lib_v2.1.0_230511.RData')
# usethis::use_data(ms2_dodd_lib, overwrite = TRUE)
#
################################################################################
# # ms2_dodd_lib -----------------------------------------------------------------
# load('~/Project/04_package/00_Database/DoddLib/02_ms2_lib/v230612_v2.2.0/ms2_dodd_lib_v2.2.0_230612.RData')
# usethis::use_data(ms2_dodd_lib, overwrite = TRUE)

################################################################################
# # ms2 peptide ----------------------------------------------------------------
# load('~/Project/04_package/00_Database/FiehnPeptide/fiehn_peptide_lib_230727.RData')
# usethis::use_data(fiehn_peptide_lib, overwrite = TRUE)

################################################################################
# # ms2_dodd_lib -----------------------------------------------------------------
# load('~/Project/04_package/00_Database/DoddLib/02_ms2_lib/v231206_v2.3.0/ms2_dodd_lib_v2.3.0_231206.RData')
# usethis::use_data(ms2_dodd_lib, overwrite = TRUE)

# # cpd_dodd_lib -----------------------------------------------------------------
# load('~/Project/04_package/00_Database/DoddLib/01_compound_RT_lib/04_version_231206/cpd_dodd_lib_v2.0.0_231206.RData')
# usethis::use_data(cpd_dodd_lib, overwrite = TRUE)


################################################################################
# gnps_lib ---------------------------------------------------------------------
# load('~/Project/04_package/00_Database/GNPS_ReverseLib/00_converted_database/ms2_gnps_bile_acid_lib.RData')
# gnps_bile_acid_lib <- ms2_gnps_bile_acid_lib
# usethis::use_data(gnps_bile_acid_lib)
#
# load('~/Project/04_package/00_Database/GNPS_ReverseLib/00_converted_database/ms2_gnps_acyl_amides_lib.RData')
# gnps_acyl_amides_lib <- ms2_gnps_acyl_amides_lib
# usethis::use_data(gnps_acyl_amides_lib)
#
# load('~/Project/04_package/00_Database/GNPS_ReverseLib/00_converted_database/ms2_gnps_acyl_esters_lib.RData')
# gnps_acyl_esters_lib <- ms2_gnps_acyl_esters_lib
# usethis::use_data(gnps_acyl_esters_lib)

################################################################################
# pathdb_enteropath ------------------------------------------------------------
# load('~/Project/04_package/00_Database/Enteropathway/pathdb_enteropathway_v1_240222.RData')
# usethis::use_data(pathdb_enteropathway)


################################################################################
# DoddLab ----------------------------------------------------------------------

# # v2.3.0
# load('~/Project/04_package/00_Database/DoddLib/02_ms2_lib/v231206_v2.3.0/ms2_dodd_lib_v2.3.0_231206.RData')
# list_ms2_dodd_lib <- list('v2.3.0' = ms2_dodd_lib)
#
# # v2.4.0
# load('~/Project/04_package/00_Database/DoddLib/02_ms2_lib/v240325_v2.4.0/ms2_dodd_lib_v2.4.0_240325.RData')
# usethis::use_data(ms2_dodd_lib, overwrite = TRUE)
#
# list_ms2_dodd_lib$v2.4.0 <- ms2_dodd_lib
# usethis::use_data(list_ms2_dodd_lib, overwrite = TRUE)


# # cpd_dodd_lib -----------------------------------------------------------------
# #
# load('~/Project/04_package/00_Database/DoddLib/01_compound_RT_lib/04_version_231206/cpd_dodd_lib_v2.0.0_231206.RData')
# list_cpd_dodd_lib <- list('v2.0.0' = cpd_dodd_lib)
#
# load('~/Project/04_package/00_Database/DoddLib/01_compound_RT_lib/08_version_240325/cpd_dodd_lib_v2.4.0_240325.RData')
# usethis::use_data(cpd_dodd_lib, overwrite = TRUE)
#
# list_cpd_dodd_lib$v2.4.0 <- cpd_dodd_lib
#
# usethis::use_data(list_cpd_dodd_lib, overwrite = TRUE)


################################################################################
# DoddLab ----------------------------------------------------------------------

# # v2.5.0
# load('~/Project/04_package/00_Database/DoddLib/02_ms2_lib/v240327_v2.5.0/ms2_dodd_lib_v2.5.0_240327.RData')
# usethis::use_data(ms2_dodd_lib, overwrite = TRUE)
#
# load('~/Project/04_package/DoddLabDatabase/data/list_ms2_dodd_lib.rda')
# save(list_ms2_dodd_lib,
#      file = '~/Project/04_package/00_Database/DoddLib/02_ms2_lib/v240325_v2.4.0/list_ms2_dodd_lib_v2.4.0_240325.RData')
#
# list_ms2_dodd_lib$v2.5.0 <- ms2_dodd_lib
# usethis::use_data(list_ms2_dodd_lib, overwrite = TRUE)
# save(list_ms2_dodd_lib,
#      file = '~/Project/04_package/00_Database/DoddLib/02_ms2_lib/v240327_v2.5.0/list_ms2_dodd_lib_v2.5.0_240327.RData')

# cpd_dodd_lib -----------------------------------------------------------------
#
#
# load('~/Project/04_package/00_Database/DoddLib/01_compound_RT_lib/10_version_240329/cpd_dodd_lib_v2.5.0_240329.RData')
# usethis::use_data(cpd_dodd_lib, overwrite = TRUE)
#
# load('~/Project/04_package/DoddLabDatabase/data/list_cpd_dodd_lib.rda')
# save(list_cpd_dodd_lib,
#      file = '~/Project/04_package/00_Database/DoddLib/01_compound_RT_lib/08_version_240325/list_cpd_dodd_lib_v2.4.0_240325.RData')
#
#
# list_cpd_dodd_lib$v2.5.0 <- cpd_dodd_lib
# usethis::use_data(list_cpd_dodd_lib, overwrite = TRUE)
#
# save(list_cpd_dodd_lib,
#      file = '~/Project/04_package/00_Database/DoddLib/01_compound_RT_lib/10_version_240329/list_cpd_dodd_lib_v2.5.0_240329.RData')

################################################################################
# Mass-DIAL lipid blast lib ----------------------------------------------------
#
# load('~/Project/04_package/00_Database/MSDIAL_lipidBlast/msdial_lipidblast_lib_240424.RData')
# usethis::use_data(msdial_lipidblast_lib, overwrite = TRUE)

################################################################################
# DoddLab ----------------------------------------------------------------------
#
# load('~/Project/04_package/00_Database/DoddLib/02_ms2_lib/v240612_v2.6.1/ms2_dodd_lib_v2.6.1_240612.RData')
# ms2_dodd_lib_v2.6.1 <- ms2_dodd_lib
# load('./data/list_ms2_dodd_lib.rda')
# list_ms2_dodd_lib$v2.6.1 <- ms2_dodd_lib_v2.6.1
# usethis::use_data(list_ms2_dodd_lib, overwrite = TRUE)
#
# save(list_ms2_dodd_lib,
#      file = '~/Project/04_package/00_Database/DoddLib/02_ms2_lib/v240612_v2.6.1/list_ms2_dodd_lib_v2.6.1_240625.RData')
#
# load('~/Project/04_package/00_Database/DoddLib/01_compound_RT_lib/12_version_240531/cpd_dodd_lib_v2.6.1_240612.RData')
# load('~/Project/04_package/00_Database/DoddLib/01_compound_RT_lib/10_version_240329/list_cpd_dodd_lib_v2.5.0_240329.RData')
#
# list_cpd_dodd_lib$v2.6.1 <- cpd_dodd_lib
# usethis::use_data(list_cpd_dodd_lib, overwrite = TRUE)
# save(list_cpd_dodd_lib,
#      file = '~/Project/04_package/00_Database/DoddLib/01_compound_RT_lib/12_version_240531/list_cpd_dodd_lib_v2.6.1_240625.RData')

################################################################################
