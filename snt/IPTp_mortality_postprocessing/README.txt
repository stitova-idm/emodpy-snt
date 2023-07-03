After DTK simulations have completed, still need to calculate the effects of MiP, IPTp, IPTi, and mortality in post-processing. 

After running simAdjustments_mortality_MiP_IPTp_IPTi.R, an output csv file is generated in each scenario's folder called malariaBurden_withAdjustments.csv.
Description of select columns in output file (stars next to those relevant for plotting output):
  *Statistical_Population - among all ages (original simulation output)
  PfHRP2_Prevalence - original simulation output, among all ages, does not include adjustments for IPTp
  Received_Severe_treatment - original simulation output - number of individuals with severe disease who were treated, does not include treatment for severe disease from MiP
  *New_Clinical_Cases - adjusted for IPTi (no change from MiP or IPTp)
  New_Severe_Cases - original simulation output - does not include severe disease from MiP
  *PfPR_MiP_adjusted - PfPR for all ages, adjusted to account for IPTp (and IPTi if relevant)
  severe_maternal - expected number of individuals with a new occurrence of severe anemia due to MiP
  severe_total - sum of 'New Severe Cases' from original simulation output and severe_maternal. represents number of individuals with severe disease due to malaria for the entire population.
  severe_maternal_treated - expected number of individuals treated among those with a new occurrence of severe anemia due to MiP
  mLBW_births - number of infants expected to be born with low birth rate due to MiP
  mLBW_deaths - number of mLBW infants expected to die within their first year who would have survived in the absence of MiP
  *MiP_stillbirths - number of stillbirths attributable to MiP with assumed IPTp coverage
  MiP_stillbirths_noIPTp - number of stillbirths attributable to MiP if no IPTp were used in the population
  direct_mortality_nonMiP_1 - number of deaths due to severe disease using larger estimate for CFR for untreated severe cases (includes mortality from severe maternal anemia, does not include mLBW or stillbirths)
  direct_mortality_nonMiP_2 - number of deaths due to severe disease using smaller estimate for CFR for untreated severe cases (includes mortality from severe maternal anemia, does not include mLBW or stillbirths)
  indirect_mortality_nonMiP_1 - number of deaths that are not directly caused by malaria but that would not have occurred in the absence of the malaria infection (e.g., due to coinfection, malnutrition, etc. does not include mLBW). uses larger estimate for CFR
  indirect_mortality_nonMiP_2 - number of deaths that are not directly caused by malaria but that would not have occurred in the absence of the malaria infection (e.g., due to coinfection, malnutrition, etc. does not include mLBW). uses smaller estimate for CFR
  *total_mortality_1 - number of deaths in population directly or indirectly attributable to malaria (includes severe materal anemia, deaths due to mLBW, direct, and indirect malaria deaths. does not include stillbirth). uses larger estimates for CFRs.
  *total_mortality_2 - number of deaths in population directly or indirectly attributable to malaria (includes severe materal anemia, deaths due to mLBW, direct, and indirect malaria deaths. does not include stillbirth). uses smaller estimates for CFRs.
  *PfPR_U5 - PfPR among U5, adjusted for IPTi
  *Pop_U5 - population size among individuals U5 (original simulation output)
  Num_U5_Received_Severe_Treatment - number of U5 individuals with severe disease who were treated from simulation (does not include mLBW or stillbirths)
  *New_clinical_cases_U5 - number of U5 individuals with clinical cases (does not include mLBW or stillbirths), adjusted for IPTi
  Severe_cases_U5 - number of U5 individuals with severe disease (does not include mLBW or stillbirths)
  direct_mortality_nonMiP_U5_1 - number of deaths in U5 due to severe disease using larger estimate for CFR for untreated severe cases (does not include mLBW or stillbirths)
  direct_mortality_nonMiP_U5_2 - number of deaths in U5 due to severe disease using smaller estimate for CFR for untreated severe cases (does not include mLBW or stillbirths)
  indirect_mortality_nonMiP_U5_1 - number of deaths in U5 that are not directly caused by malaria but that would not have occurred in the absence of the malaria infection (e.g., due to coinfection, malnutrition, etc. does not include mLBW). uses larger estimate for CFR
  indirect_mortality_nonMiP_U5_2 - number of deaths in U5 that are not directly caused by malaria but that would not have occurred in the absence of the malaria infection (e.g., due to coinfection, malnutrition, etc. does not include mLBW). uses smaller estimate for CFR
  *total_mortality_U5_1 - number of deaths in U5 directly or indirectly attributable to malaria (includes deaths due to mLBW, direct, and indirect malaria deaths. does not include stillbirth). uses larger estimates for CFRs.
  *total_mortality_U5_2 - number of deaths in U5 directly or indirectly attributable to malaria (includes deaths due to mLBW, direct, and indirect malaria deaths. does not include stillbirth). uses smaller estimates for CFRs.
