
********************************************
// fairness_sem: Multiple Group Test of 
//		Schema Invariance for Fairness paper
// Note: Forthcomin in SER 
// Author: Marshall A. Taylor
********************************************

version 16.0
clear all
macro drop _all
log using fairness_sem.log, replace text
set more off  


**************************
***COMMANDS BEGIN HERE ***
**************************

//Load in data
import delimited "data/sem_data.csv", clear


//Encode grouping variable and remove NA category
encode group, gen(group2)
recode group2 (5 = .), gen(group3)


//SEMs
sem painter autoa ///
			redwageb fixedinc pb grocery ///
			apples dolls hiring, group(group3)
est store DIFF9

sem painter autoa ///
			redwageb fixedinc pb grocery ///
			apples dolls hiring, group(group3) ginvariant(covex)
est store SAME9

est stats DIFF9 SAME9
lrtest SAME9 DIFF9

**************************
*** COMMANDS END HERE  ***
**************************
log close
exit
