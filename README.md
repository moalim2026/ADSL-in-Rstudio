# ADSL-in-Rstudio
Creating ADSL in R

library(admiral)
library(dplyr, warn.conflicts = FALSE)
library(pharmaversesdtm)
library(lubridate)
library(stringr)

library(tibble) # added 


dm <- pharmaversesdtm::dm
ds <- pharmaversesdtm::ds
ex <- pharmaversesdtm::ex
ae <- pharmaversesdtm::ae
lb <- pharmaversesdtm::lb

dm <- convert_blanks_to_na(dm)
ds <- convert_blanks_to_na(ds)
ex <- convert_blanks_to_na(ex)
ae <- convert_blanks_to_na(ae)
lb <- convert_blanks_to_na(lb)


#Derive Treatment Variables (TRT0xP, TRT0xA)
adsl <- dm %>%
  mutate(TRT01P = ARM, TRT01A = ACTARM)

#Derive/Impute Numeric Treatment Date/Time and Duration (TRTSDTM, TRTEDTM, TRTDURD)

# Impute start and end time of exposure to first and last respectively,
# Do not impute date
#Conversion and imputation is done by derive_vars_dtm().
ex_ext <- ex %>%
  derive_vars_dtm(
    dtc = EXSTDTC,
    new_vars_prefix = "EXST"
  ) %>%
  derive_vars_dtm(
    dtc = EXENDTC,
    new_vars_prefix = "EXEN",
    time_imputation = "last"
  )

adsl <- adsl %>%
  derive_vars_merged(
    dataset_add = ex_ext,
    filter_add = (EXDOSE > 0 |
                    (EXDOSE == 0 &
                       str_detect(EXTRT, "PLACEBO"))) & !is.na(EXSTDTM),
    new_vars = exprs(TRTSDTM = EXSTDTM, TRTSTMF = EXSTTMF),
    order = exprs(EXSTDTM, EXSEQ),
    mode = "first",
    by_vars = exprs(STUDYID, USUBJID)
  ) %>%
  derive_vars_merged(
    dataset_add = ex_ext,
    filter_add = (EXDOSE > 0 |
                    (EXDOSE == 0 &
                       str_detect(EXTRT, "PLACEBO"))) & !is.na(EXENDTM),
    new_vars = exprs(TRTEDTM = EXENDTM, TRTETMF = EXENTMF),
    order = exprs(EXENDTM, EXSEQ),
    mode = "last",
    by_vars = exprs(STUDYID, USUBJID)
  )

#The datetime variables returned can be converted to dates using the derive_vars_dtm_to_dt() function.
adsl <- adsl %>%
  derive_vars_dtm_to_dt(source_vars = exprs(TRTSDTM, TRTEDTM))

#Now, that TRTSDT and TRTEDT are derived, the function derive_var_trtdurd() 
#can be used to calculate the Treatment duration (TRTDURD).

adsl <- adsl %>%
  derive_var_trtdurd()


# Derive Disposition Variables
# Disposition Dates (e.g. EOSDT)
# Convert character date to numeric date without imputation
# using derive_vars_dt() function
ds_ext <- derive_vars_dt(
  ds,
  dtc = DSSTDTC,
  new_vars_prefix = "DSST"
)

#To add the End of Study date (EOSDT) to the input dataset, a call could be:

adsl <- adsl %>%
  derive_vars_merged(
    dataset_add = ds_ext,
    by_vars = exprs(STUDYID, USUBJID),
    new_vars = exprs(EOSDT = DSSTDT),
    filter_add = DSCAT == "DISPOSITION EVENT" & DSDECOD != "SCREEN FAILURE"
  )

# Disposition Status (e.g. EOSSTT)
# formatting eosstt
format_eosstt <- function(x) {
  case_when(
    x %in% c("COMPLETED") ~ "COMPLETED",
    x %in% c("SCREEN FAILURE") ~ NA_character_,
    TRUE ~ "DISCONTINUED"
  )
}

#The customized mapping function format_eosstt() can now be passed to the main function. 
#For subjects without a disposition event the end of study status is set to "ONGOING" 
#by specifying the missing_values argument.

adsl <- adsl %>%  #This call would return the input dataset with the column EOSSTT added.
  
  derive_vars_merged(
    dataset_add = ds,
    by_vars = exprs(STUDYID, USUBJID),
    filter_add = DSCAT == "DISPOSITION EVENT",
    new_vars = exprs(EOSSTT = format_eosstt(DSDECOD)),
    missing_values = exprs(EOSSTT = "ONGOING")
  )

#Disposition Reason(s) (e.g. DCSREAS, DCSREASP)
#To derive the End of Study reason(s) (DCSREAS and DCSREASP), the function will map DCSREAS as DSDECOD, 
#and DCSREASP as DSTERM if DSDECOD is not "COMPLETED", "SCREEN FAILURE", or NA, NA otherwise.

adsl <- adsl %>%
  derive_vars_merged(
    dataset_add = ds,
    by_vars = exprs(USUBJID),
    new_vars = exprs(DCSREAS = DSDECOD),
    filter_add = DSCAT == "DISPOSITION EVENT" &
      DSDECOD %notin% c("SCREEN FAILURE", "COMPLETED", NA)
  ) %>%
  derive_vars_merged(
    dataset_add = ds,
    by_vars = exprs(USUBJID),
    new_vars = exprs(DCSREASP = DSTERM),
    filter_add = DSCAT == "DISPOSITION EVENT" & DSDECOD %in% "OTHER"
  )

#Randomization Date (RANDDT)
#The function derive_vars_merged() can be used to derive randomization date variable. To map Randomization Date (RANDDT), the call would be:
  
  adsl <- adsl %>%
  derive_vars_merged(
    dataset_add = ds_ext,
    filter_add = DSDECOD == "RANDOMIZED",
    by_vars = exprs(STUDYID, USUBJID),
    new_vars = exprs(RANDDT = DSSTDT)
  )
  
#Derive Birth Date and Analysis Age (BRTHDT, AAGE, AAGEU)
#Note that BRTHDT must be derived first from BRTHDTC using derive_vars_dt().
  
  # Derive birth date from BRTHDTC
  adsl <- adsl %>%
    derive_vars_dt(
      new_vars_prefix = "BRTH",
      dtc = BRTHDTC
    )

# Now that we have BRTHDT, we can use derive_var_aage() to derive AAGE and AAGEU.
  #This call returns the input dataset with AAGE and AAGEU added. By default, 
  #the age is calculated in years.
  
  adsl <- adsl %>%
    derive_vars_aage(
      start_date = BRTHDT,
      end_date = RANDDT
    )
#Cause of Death (DTHCAUS)
#For example, if the date of death is collected in the AE form when the AE is Fatal, 
#the cause of death would be set to the preferred term (AEDECOD) of that Fatal AE, 
  #while if the date of death is collected in the DS form, the cause of death would be set 
  #to the disposition term (DSTERM). To achieve this, the event() 
  #objects within derive_vars_extreme_event() must be specified and defined such that they fit 
  #the study requirement.
  
  #An example call to derive_vars_extreme_event() would be:
  adsl <- adsl %>%
    derive_vars_extreme_event(
      by_vars = exprs(STUDYID, USUBJID),
      events = list(
        event(
          dataset_name = "ae",
          condition = AEOUT == "FATAL",
          set_values_to = exprs(DTHCAUS = AEDECOD),
        ),
        event(
          dataset_name = "ds",
          condition = DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM),
          set_values_to = exprs(DTHCAUS = DSTERM),
        )
      ),
      source_datasets = list(ae = ae, ds = ds),
      tmp_event_nr_var = event_nr,
      order = exprs(event_nr),
      mode = "first",
      new_vars = exprs(DTHCAUS)
    )
  
#The function also offers the option to add some traceability variables (e.g. DTHDOM would store the 
  #domain where the date of death is collected, and DTHSEQ would store the xxSEQ value of that domain).
  #The traceability variables should be added to the event() calls and included in the new_vars 
  #parameter of derive_vars_extreme_event().
  
  
  adsl <- adsl %>%
    select(-DTHCAUS) %>% # Remove it before deriving it again
    derive_vars_extreme_event(
      by_vars = exprs(STUDYID, USUBJID),
      events = list(
        event(
          dataset_name = "ae",
          condition = AEOUT == "FATAL",
          set_values_to = exprs(DTHCAUS = AEDECOD, DTHDOM = "AE", DTHSEQ = AESEQ),
        ),
        event(
          dataset_name = "ds",
          condition = DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM),
          set_values_to = exprs(DTHCAUS = DSTERM, DTHDOM = "DS", DTHSEQ = DSSEQ),
        )
      ),
      source_datasets = list(ae = ae, ds = ds),
      tmp_event_nr_var = event_nr,
      order = exprs(event_nr),
      mode = "first",
      new_vars = exprs(DTHCAUS, DTHDOM, DTHSEQ)
    )
  
  
  
  #Derive the Last Date Known Alive (LSTALVDT)
  #Similarly as for the cause of death (DTHCAUS), the last known alive date (LSTALVDT) can be derived 
  #from multiples sources using derive_vars_extreme_event().An example could be (--DTC dates are
  #converted to numeric dates imputing missing day and month to the first):
  
  adsl <- adsl %>%
    derive_vars_extreme_event(
      by_vars = exprs(STUDYID, USUBJID),
      events = list(
        event(
          dataset_name = "ae",
          order = exprs(AESTDTC, AESEQ),
          condition = !is.na(AESTDTC),
          set_values_to = exprs(
            LSTALVDT = convert_dtc_to_dt(AESTDTC, highest_imputation = "M"),
            seq = AESEQ
          ),
        ),
        event(
          dataset_name = "ae",
          order = exprs(AEENDTC, AESEQ),
          condition = !is.na(AEENDTC),
          set_values_to = exprs(
            LSTALVDT = convert_dtc_to_dt(AEENDTC, highest_imputation = "M"),
            seq = AESEQ
          ),
        ),
        event(
          dataset_name = "lb",
          order = exprs(LBDTC, LBSEQ),
          condition = !is.na(LBDTC),
          set_values_to = exprs(
            LSTALVDT = convert_dtc_to_dt(LBDTC, highest_imputation = "M"),
            seq = LBSEQ
          ),
        ),
        event(
          dataset_name = "adsl",
          condition = !is.na(TRTEDT),
          set_values_to = exprs(LSTALVDT = TRTEDT, seq = 0),
        )
      ),
      source_datasets = list(ae = ae, lb = lb, adsl = adsl),
      tmp_event_nr_var = event_nr,
      order = exprs(LSTALVDT, seq, event_nr),
      mode = "last",
      new_vars = exprs(LSTALVDT)
    )
  
#Traceability variables can be added by specifying the variables in the set_values_to parameter 
  #of the event() function.
  
  adsl <- adsl %>%
    select(-LSTALVDT) %>% # Created in the previous call
    derive_vars_extreme_event(
      by_vars = exprs(STUDYID, USUBJID),
      events = list(
        event(
          dataset_name = "ae",
          order = exprs(AESTDTC, AESEQ),
          condition = !is.na(AESTDTC),
          set_values_to = exprs(
            LSTALVDT = convert_dtc_to_dt(AESTDTC, highest_imputation = "M"),
            LALVSEQ = AESEQ,
            LALVDOM = "AE",
            LALVVAR = "AESTDTC"
          ),
        ),
        event(
          dataset_name = "ae",
          order = exprs(AEENDTC, AESEQ),
          condition = !is.na(AEENDTC),
          set_values_to = exprs(
            LSTALVDT = convert_dtc_to_dt(AEENDTC, highest_imputation = "M"),
            LALVSEQ = AESEQ,
            LALVDOM = "AE",
            LALVVAR = "AEENDTC"
          ),
        ),
        event(
          dataset_name = "lb",
          order = exprs(LBDTC, LBSEQ),
          condition = !is.na(LBDTC),
          set_values_to = exprs(
            LSTALVDT = convert_dtc_to_dt(LBDTC, highest_imputation = "M"),
            LALVSEQ = LBSEQ,
            LALVDOM = "LB",
            LALVVAR = "LBDTC"
          ),
        ),
        event(
          dataset_name = "adsl",
          condition = !is.na(TRTEDT),
          set_values_to = exprs(LSTALVDT = TRTEDT, LALVSEQ = NA_integer_, LALVDOM = "ADSL", LALVVAR = "TRTEDTM"),
        )
      ),
      source_datasets = list(ae = ae, lb = lb, adsl = adsl),
      tmp_event_nr_var = event_nr,
      order = exprs(LSTALVDT, LALVSEQ, event_nr),
      mode = "last",
      new_vars = exprs(LSTALVDT, LALVSEQ, LALVDOM, LALVVAR)
    )
  
#Derive Groupings and Populations
#Grouping (e.g. AGEGR1 or REGION1)
  # Create lookup tables
  agegr1_lookup <- exprs(
    ~condition,           ~AGEGR1,
    AGE < 18,               "<18",
    between(AGE, 18, 64), "18-64",
    AGE > 64,               ">64",
    is.na(AGE),         "Missing"
  )
  
  region1_lookup <- exprs(
    ~condition,                          ~REGION1,
    COUNTRY %in% c("CAN", "USA"), "North America",
    !is.na(COUNTRY),          "Rest of the World",
    is.na(COUNTRY),                     "Missing"
  )
  adsl <- adsl %>%
    derive_vars_cat(
      definition = agegr1_lookup
    ) %>%
    derive_vars_cat(
      definition = region1_lookup
    )
#Alternatively, you can also solve this task with custom functions:
  
  format_agegr1 <- function(var_input) {
    case_when(
      var_input < 18 ~ "<18",
      between(var_input, 18, 64) ~ "18-64",
      var_input > 64 ~ ">64",
      TRUE ~ "Missing"
    )
  }
  format_region1 <- function(var_input) {
    case_when(
      var_input %in% c("CAN", "USA") ~ "North America",
      !is.na(var_input) ~ "Rest of the World",
      TRUE ~ "Missing"
    )
  }
  
  adsl %>%
    mutate(
      AGEGR1 = format_agegr1(AAGE),
      REGION1 = format_region1(COUNTRY)
    )
  
  
  #Population Flags (e.g. SAFFL)
  #Since the populations flags are mainly company/study specific no dedicated functions are provided, 
  #but in most cases they can easily be derived using derive_var_merged_exist_flag.
  
  #An example of an implementation could be:
    
    adsl <- adsl %>%
    derive_var_merged_exist_flag(
      dataset_add = ex,
      by_vars = exprs(STUDYID, USUBJID),
      new_var = SAFFL,
      false_value = "N",
      missing_value = "N",
      condition = (EXDOSE > 0 | (EXDOSE == 0 & str_detect(EXTRT, "PLACEBO")))
    )
    
    
    # Select specific columns using base R
    selected_columns <- adsl[, c("USUBJID","AGEGR1", "REGION1", "TRT01P")]
    print(selected_columns)

   

    
    
    
    
    
