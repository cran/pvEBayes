#' FDA statin dataset with 44 adverse events
#'
#' An adverse event-drug count dataset (contingency table)
#' obtained from the FDA FAERS database for the
#' quarters 2021Q1 - 2024Q4.
#'
#'
#' @details
#' Data are stored in the form of a contingency table, with
#' drugs listed across the columns and the 44 AEs presented across
#' the rows. Each cell in the contingency table represents the total
#' report counts associated with that (AE, drug) pair and detected
#' in the FDA FAERS database during 2021Q1 - 2024Q4.
#'
#' The dataset catalogs 6 statin drugs (across columns):
#'
#' Atorvastatin, Fluvastatin, Lovastatin, Pravastatin, Rosuvastatin, Simvastatin.
#'
#' The 44 adverse events presented across the rows are considered
#' significant by FDA.
#'
#' This dataset is an updated version of statin46 from the pvLRT package which
#' collect the same scope of AEs for 6 statin drugs for quarters 2014Q3 - 2020Q4.
#'
#' During 2021Q1 - 2024Q4, there was no AE report for "BLOOD CREATINE
#' PHOSPHOKINASE MM INCREASED" and "MYOGLOBIN BLOOD PRESENT". Therefore, these
#' two AEs are not presented in the statin2025_44 dataset.
#'
#'
#'
#'
#'
#'
#' @source \url{https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html}
"statin2025_44"



#' FDA statin dataset with 42 adverse events
#'
#' An adverse event-drug count dataset (contingency table)
#' obtained from the FDA FAERS database for the
#' quarters 2014Q3 - 2020Q4.
#'
#'
#' @details
#' Data are stored in the form of a contingency table, with
#' drugs listed across the columns and the 42 AEs presented across
#' the rows. Each cell in the contingency table represents the total
#' report counts associated with that (AE, drug) pair and detected
#' in the FDA FAERS database during 2014Q3 - 2020Q4.
#'
#' The dataset catalogs 6 statin drugs (across columns):
#'
#' Atorvastatin, Fluvastatin, Lovastatin, Pravastatin, Rosuvastatin, Simvastatin.
#'
#'
#' This dataset is derived from the `statin46` dataset in the \pkg{pvLRT}
#' package, with four AEs removed.
#'
#' The excluded AEs are:
#' "Blood Creatine Phosphokinase Mm Increased",
#' "Myoglobin Blood Present",
#' "Myoglobin Urine Present", and
#' "Myoglobinaemia".
#'
#'
#'
#'
#'
#'
#'
#' @source \url{https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html}
"statin42"











#' FDA statin dataset with 5119 adverse events
#'
#' An adverse event-drug count dataset (contingency table)
#' obtained from the FDA FAERS database for the
#' quarters 2021Q1 - 2024Q4.
#'
#'
#' @details
#'
#' The dataset catalogs 6 statin drugs (across columns):
#' Atorvastatin, Fluvastatin, Lovastatin, Pravastatin, Rosuvastatin, Simvastatin.
#'
#' Data are stored in the form of a contingency table, with
#' drugs listed across the columns and the 5119 AEs presented across
#' the rows. Each cell in the contingency table represents the total
#' report counts associated with that (AE, drug) pair and detected
#' in the FDA FAERS database during 2021Q1 - 2024Q4.
#'
#' The dataset catalogs 6 statin drugs (across columns):
#'
#' Atorvastatin, Fluvastatin, Lovastatin, Pravastatin, Rosuvastatin, Simvastatin.
#'
#' The 5119 adverse events presented across the rows are AEs that contain at
#' least one report for 6 statin drugs during 2021Q1 - 2024Q4.
#'
#' This dataset is an updated version of statin from the pvLRT package which collects
#' the same scope of AEs for 6 statin drugs for quarters 2014Q3 - 2020Q4.
#'
#'
#'
#'
#' @source \url{https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html}
"statin2025"




#' FDA GBCA dataset with 1328 adverse events
#'
#' An adverse event-drug count dataset (contingency table)
#' obtained from the FDA FAERS database for the
#' quarters 2021Q1 - 2024Q4.
#'
#'
#' @details
#' Data are stored in the form of a contingency table, with
#' drugs listed across the columns and the 1328 AEs presented across
#' the rows. Each cell in the contingency table represents the total
#' report counts associated with that (AE, drug) pair and detected
#' in the FDA FAERS database during 2021Q1 - 2024Q4.
#'
#' The dataset catalogs 7 Gadolinium-Based Contrast Agents (GBCAs) (across columns):
#'
#' Gadobenate, Gadobutrol, Gadodiamide, Gadopentetate, Gadoterate, Gadoteridol,
#' Gadoxetate
#'
#' The 1328 adverse events presented across the rows are AEs that contain at
#' least one report for the 7 GBCA drugs during 2021Q1 - 2024Q4.
#'
#' This dataset is an updated version of gbca from the pvLRT package which collects
#' the same scope of AEs for 7 gbca drugs for quarters 2014Q3 - 2020Q4.
#'
#'
#'
#'
#' @source \url{https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html}
"gbca2025"



#' FDA GBCA dataset with 69 adverse events
#'
#' An adverse event-drug count dataset (contingency table)
#' obtained from the FDA FAERS database for the
#' quarters 2021Q1 - 2024Q4
#'
#'
#' @details
#' Data are stored in the form of a contingency table, with
#' drugs listed across the columns and the 69 AEs presented across
#' the rows. Each cell in the contingency table represents the total
#' report counts associated with that (AE, drug) pair and detected
#' in the FDA FAERS database during 2021Q1 - 2024Q4.
#'
#' The dataset catalogs 7 Gadolinium-Based Contrast Agents (GBCAs)
#' (across columns):
#'
#' Gadobenate, Gadobutrol, Gadodiamide, Gadopentetate, Gadoterate, Gadoteridol,
#' Gadoxetate.
#'
#' The 69 adverse events presented across the rows are selected from 1328 AEs of
#' gbca2025 which are related to the brain or neural system. Other AEs are collapsed
#' to one reference row: "Other AEs".
#'
#'
#'
#'
#'
#'
#' @source \url{https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html}
"gbca2025_69"
