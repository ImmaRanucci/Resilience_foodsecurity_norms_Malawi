/*

Title: Weather Shocks and Resilience to Food Insecurity: Exploring the Role of Gender and Kinship Norms

Authors: Immacolata Ranucci; Donato Romano; Luca Tiberti

Replication code


*******************************************************************************/

*** Filepath 

clear all
	gl  PATH "C:\Users\USERNAME\git\Resilience_foodsecurity_norms_Malawi"
	cd "${PATH}"


capture mkdir "${PATH}\Results"

capture mkdir "${PATH}\Results\Tables"
gl TAB ${PATH}\Results\Tables

capture mkdir "${PATH}\Results\Figures"
gl FIG ${PATH}\Results\Figures

*_______________________________________________________________________________

*** Install packages

ssc install estout, all replace
ssc install grc1leg2, replace
ssc install exbsample, replace

*_______________________________________________________________________________


***	Load data
use "$PATH/Data.dta", clear


*_______________________________________________________________________________

*								SAMPLE SELECTION 

* Table 1: Number of observations after the application of each selection criteria
tab year				// 0. Original panel 
tab year if sample1==1  // 1. Rural households 
tab year if sample2==1  // 2. Agricultural households
tab year if sample3==1  // 3. Consumption recorded by 31Octb
tab year if sample4==1  // 4. Observed at least in two consecutive rounds
tab year if sample5==1  // 5. No Joint Land Managementd


*_______________________________________________________________________________

*							SUMMARY STATISTICS

		
keep if sample4==1

* Table 2: Summary statistics
	eststo clear
	macro def SUM FCS HDDS spei_6 drought flood spei_6_lag droughtlag lm1fem lm1joint lm12fem lm12joint comm_marriagetype_matri hhsize dependency  fhead age_head edu_adult split1013 split1316 incmlft hh_diffcomm2010 wealth_index arable_land cropdiv_SDI intercrop TLU ssn_freemaize ssn_freefood ssn_FCfW ssn_childfeed
	eststo: quietly estpost summarize $SUM [weight=panelwgt] if year==2010
	eststo: quietly estpost summarize $SUM [weight=panelwgt] if year==2013 
	eststo: quietly estpost summarize $SUM [weight=panelwgt] if year==2016
	eststo: quietly estpost summarize $SUM [weight=panelwgt] 
	esttab using  "$TAB/Table2_sumstat.tex", tex cells("count(fmt(0))  mean(fmt(2)) sd(fmt(2))") label nodepvar replace title("Table 2: Summary statistics")


*_______________________________________________________________________________

*** 							ANALYSIS

keep if sample5==1

*** Panel setting
	xtset HHID y
	xtdes
	
*** Create replication weights for bootstrap - Exponential distribution
	exbsample 500 [pweight=panelwgt] , stub(bwexp) dis(exponential) cluster(EA parent_hh) strata(stratum) seed(220796)

*** Survey setting with bootstrap weights
	svyset EA , strata(stratum)  || HHID [pw=panelwgt], bsrweight(bw*)


*_______________________________________________________________________________

*						Estimation of Equation 1

	gen W = .
	gen Wsq=.
	gen Wtr=.
	gen thresh=.
	
*** Definition of macros for variables
	macro def MM 	i.comm_marriagetype_matri 
	macro def INT 	i.drought##i.landmanagement1
	macro def T 	i.y i.consumption_month i.split1013 i.split1316		
	macro def G 	i.grid
	global INDEP 	$INT $T $G
	

*-------------------------------------------------------------------------------
* 						

*** FOOD CONSUMPTION SCORE 

	* Outcome
		replace W=FCS
		replace Wsq = W*W
		replace Wtr = W*W*W
        macro def LAG c.l.W c.l.Wsq c.l.Wtr

	* Drop variables
			capture drop resid_cons* cond_mean  cond_var* prob* alfa beta  
			svy bootstrap , nodrop : reg W ($LAG $INDEP )##$MM 
				estimate store FCS
				gen sample_analysis=1 if e(sample)
				sum FCS
				estadd scalar Wmean = r(mean): FCS
		* Conditional mean: 
				predict cond_mean, xb 
		* Conditional variance: 
				predict resid_cons, residuals
				gen resid_cons_sq_log = ln(resid_cons*resid_cons +1)
			svy bootstrap , nodrop : reg resid_cons_sq_log ($LAG $INDEP )##$MM 
				predict cond_var_log, xb
				gen cond_var= exp(cond_var_log) - 1
		* Setting Parameters for Gamma distribution
				gen alfa=(cond_mean*cond_mean)/cond_var
				gen beta=cond_var/cond_mean
		* Resilience scores:
			foreach num of numlist 35 44 {
			replace thresh= `num'
			gen prob_above_`num'=1- gammap(alfa, thresh/beta)
			svy bootstrap, nodrop : reg prob_above_`num' ($LAG $INDEP )##$MM 
				estimate store FCS_RES`num'
				sum prob_above_`num'
				estadd scalar Wmean = r(mean): FCS_RES`num'
			gen RS_FCS_`num'=prob_above_`num'
			}
			
		
			
*-------------------------------------------------------------------------------

*** HOUSEHOLD DIETARY DIVERSITY SCORE - Normal distribution

	* Outcome
		replace W= HDDS
		replace Wsq = W*W
		replace Wtr = W*W*W
        macro def LAG c.l.W c.l.Wsq c.l.Wtr

	*  Drop variables
			drop resid_cons* cond_mean  cond_var*  prob*  alfa beta // 
			svy bootstrap , nodrop : reg W ($LAG $INDEP )##$MM 
				estimate store HDDS
				sum HDDS
				estadd scalar Wmean = r(mean): HDDS
		* Conditional mean: 
				predict cond_mean, xb 
		* Conditional variance: 
				predict resid_cons, residuals
				gen resid_cons_sq_log = ln(resid_cons*resid_cons +1)
			svy bootstrap , nodrop : reg resid_cons_sq_log ($LAG $INDEP )##$MM 
				predict cond_var_log, xb
				gen cond_var= exp(cond_var_log) - 1
		* Setting Parameters for Gamma distribution
				gen alfa=(cond_mean*cond_mean)/cond_var
				gen beta=cond_var/cond_mean
		* Resilience scores:
			foreach num of numlist 5 7 {
			replace thresh= `num'
			gen prob_above_`num'=1- gammap(alfa, thresh/beta)
			svy bootstrap , nodrop :  reg prob_above_`num' ($LAG $INDEP )##$MM 
				estimate store HDDS_RES`num'
				sum prob_above_`num'
				estadd scalar Wmean = r(mean): HDDS_RES`num'
			gen RS_HDDS_`num'=prob_above_`num'
			}

drop bwexp*

*-------------------------------------------------------------------------------

*** Table 3: Estimation of Equation 1

			estout  FCS  FCS_RES35 FCS_RES44 HDDS  HDDS_RES5 HDDS_RES7 , ///
					cells(b(star fmt(3)) se(par fmt(3))) starlevels( * 0.10 ** 0.05 *** 0.010) stat(N  r2 N_reps, fmt(0  3 0)) ///
					keep(1.drought 1.landmanagement1 1.drought#1.landmanagement1 1.comm_marriagetype_matri 1.comm_marriagetype_matri#cL.W  1.drought#1.comm_marriagetype_matri 1.landmanagement1#1.comm_marriagetype_matri 1.drought#1.landmanagement1#1.comm_marriagetype_matri) 
									

			esttab 	FCS  FCS_RES35 FCS_RES44 ///
					HDDS  HDDS_RES5 HDDS_RES7  ///
					using  "$TAB/Table_3_Equation1.tex",  /// 
					b(%9.3f) se(%9.3f) star(* 0.10 ** 0.05 *** 0.01) ///
					stat(Wmean N  N_reps r2 , fmt(2 0 0 3 ) ///
					label("Mean of dependent var.:" "N obs." "N. reps" "$ R^2 $"  )) ///
					tex replace nogaps nobaselevels noconstant noomitted ///
					coeflabels(1.drought "Drought" 1.landmanagement1 "Female Land Management" 1.comm_marriagetype_matri "Norm (1=Matrilineal-Matrilocal)"   1.drought#1.comm_marriagetype_matri "Drought x Norm" 1.landmanagement1#1.comm_marriagetype_matri "Female LM x Norm" 1.drought#1.landmanagement1 "Drought x Female LM" 1.drought#1.landmanagement1#1.comm_marriagetype_matri "Drought x Female LM X Norm" )  ///
					order(1.drought 1.landmanagement1 1.drought#1.landmanagement1 1.comm_marriagetype_matri 1.comm_marriagetype_matri#cL.W  1.drought#1.comm_marriagetype_matri 1.landmanagement1#1.comm_marriagetype_matri 1.drought#1.landmanagement1#1.comm_marriagetype_matri) ///
					title("Estimation results for Model \ref{m4}") ///
					mtitles( "FCS" "RS: p(FCS>35)"  "RS: p(FCS>44)" "HDDS" "RS: p(HDDS>5)"  "RS: p(HDDS>7)" ) ///
					 drop (_cons *.y* *.consumption_month* *.grid* *split* *W* ) 

*-------------------------------------------------------------------------------

*** Figure 2: Distribution of Resilience Scores by drought, land management (LM) and Norm		

qui: twoway (kdensity RS_FCS_35 if comm_marriagetype_matri==1 & landmanagement1==1 & drought==0, lc(gs0)  lp(solid) lwidth(medthick) legend(label(1 "Female LM in Matrilineal-Matrilocal communities"))) ///
       (kdensity RS_FCS_35 if comm_marriagetype_matri==0 & landmanagement1==1 & drought==0, lc(gs0)  lp(shortdash) lwidth(medthick) legend(label(2 "Female LM in Other communities"))) ///
	   (kdensity RS_FCS_35 if comm_marriagetype_matri==1 & landmanagement1==0 & drought==0, lc(gs10)  lp(solid) lwidth(medthick) legend(label(3 "Male LM in Matrilineal-Matrilocal communities")))  ///
	   (kdensity RS_FCS_35 if comm_marriagetype_matri==0 & landmanagement1==0 & drought==0, lc(gs10)  lp(shortdash) lwidth(medthick) legend(label(4 "Male LM in Other communities"))) , ///
	   legend(rows(2) size(small) symx(3pt) rowgap(0.2pt))  ysize(10) xsize(10) ylabel(, labsize(large)) xlabel(, labsize(large)) ///
	   ytitle("Density", size(large)) xtitle("RS=p(FCS>35|X)", size(large)) ysc(r(3)) ///
	   graphregion(fcolor(white)) name(RS_FCS_35_D0, replace) saving(RS_FCS_35_D0, replace) ///
	   legend(rows(4) pos(3) size(large) symx(20) rowgap(0.2pt)) legend(region( lcolor(black) )) ///
	   title(No Drought) 
qui: twoway (kdensity RS_FCS_35 if comm_marriagetype_matri==1 & landmanagement1==1 & drought==1, lc(gs0)  lp(solid) lwidth(medthick)) ///
       (kdensity RS_FCS_35 if comm_marriagetype_matri==0 & landmanagement1==1 & drought==1, lc(gs0)  lp(shortdash) lwidth(medthick)) ///
	   (kdensity RS_FCS_35 if comm_marriagetype_matri==1 & landmanagement1==0 & drought==1, lc(gs10)  lp(solid) lwidth(medthick))  ///
	   (kdensity RS_FCS_35 if comm_marriagetype_matri==0 & landmanagement1==0 & drought==1, lc(gs10)  lp(shortdash) lwidth(medthick)) , ///
	   legend(rows(2) size(small) symx(3pt) rowgap(0.2pt))  ysize(10) xsize(10) ylabel(, labsize(large)) xlabel(, labsize(large)) ///
	   ytitle("Density", size(large)) xtitle("RS=p(FCS>35|X)", size(large)) ysc(r(3)) ///
	   graphregion(fcolor(white)) name(RS_FCS_35_D1, replace) saving(RS_FCS_35_D1, replace) ///
	   title(Drought) 
grc1leg2 RS_FCS_35_D0.gph  RS_FCS_35_D1.gph , span  rows(1) xcommon ycommon graphregion(fcolor(white)) name(RSFCS35_m4, replace) saving(RSFCS35_m4, replace) xsize(10) ysize(5)  title("Food Consumption Score", size(medlarge)) subtitle("Lower threshold", size(medium))  loff

qui: twoway (kdensity RS_HDDS_7 if comm_marriagetype_matri==1 & landmanagement1==1 & drought==0, lc(gs0)  lp(solid) lwidth(medthick) legend(label(1 "Female LM in Matrilineal-Matrilocal communities"))) ///
       (kdensity RS_HDDS_7 if comm_marriagetype_matri==0 & landmanagement1==1 & drought==0, lc(gs0)  lp(shortdash) lwidth(medthick) legend(label(2 "Female LM in Other communities"))) ///
	   (kdensity RS_HDDS_7 if comm_marriagetype_matri==1 & landmanagement1==0 & drought==0, lc(gs10)  lp(solid) lwidth(medthick) legend(label(3 "Male LM in Matrilineal-Matrilocal communities")))  ///
	   (kdensity RS_HDDS_7 if comm_marriagetype_matri==0 & landmanagement1==0 & drought==0, lc(gs10)  lp(shortdash) lwidth(medthick) legend(label(4 "Male LM in Other communities"))) , ///
	   legend(rows(2) size(small) symx(3pt) rowgap(0.2pt))  ysize(10) xsize(10) ylabel(, labsize(large)) xlabel(, labsize(large)) ///
	   ytitle("Density", size(large)) xtitle("RS=p(HDDS>7|X)", size(large)) ysc(r(3)) ///
	   graphregion(fcolor(white)) name(RS_HDDS_7_D0, replace) saving(RS_HDDS_7_D0, replace) ///
	   legend(rows(4) pos(3) size(large) symx(20) rowgap(0.2pt)) legend(region( lcolor(black) )) ///
	   title(No Drought) 
qui: twoway (kdensity RS_HDDS_7 if comm_marriagetype_matri==1 & landmanagement1==1 & drought==1, lc(gs0)  lp(solid) lwidth(medthick)) ///
       (kdensity RS_HDDS_7 if comm_marriagetype_matri==0 & landmanagement1==1 & drought==1, lc(gs0)  lp(shortdash) lwidth(medthick)) ///
	   (kdensity RS_HDDS_7 if comm_marriagetype_matri==1 & landmanagement1==0 & drought==1, lc(gs10)  lp(solid) lwidth(medthick))  ///
	   (kdensity RS_HDDS_7 if comm_marriagetype_matri==0 & landmanagement1==0 & drought==1, lc(gs10)  lp(shortdash) lwidth(medthick)) , ///
	   legend(rows(2) size(small) symx(3pt) rowgap(0.2pt))  ysize(10) xsize(10) ylabel(, labsize(large)) xlabel(, labsize(large)) ///
	   ytitle("Density", size(large)) xtitle("RS=p(HDDS>7|X)", size(large)) ysc(r(3)) ///
	   graphregion(fcolor(white)) name(RS_HDDS_7_D1, replace) saving(RS_HDDS_7_D1, replace) ///
	   title(Drought) 
grc1leg2 RS_HDDS_7_D0.gph  RS_HDDS_7_D1.gph , span  rows(1) xcommon ycommon graphregion(fcolor(white)) name(RSHDDS7_m4, replace) saving(RSHDDS7_m4, replace) xsize(10) ysize(5)  title("Household Dietary Diversity Score", size(medlarge)) subtitle("Higher threshold", size(medium))   // loff

grc1leg2 RSFCS35_m4.gph  RSHDDS7_m4.gph RS_FCS_35_D0.gph , span legendfrom(RS_FCS_35_D0.gph) hidelegendfrom rows(3) graphregion(fcolor(white)) name(RS_distribution, replace) saving(RS_distribution, replace) xsize(25) ysize(30) labsize(*0.5) symxsize(*3) // loff
graph export "${FIG}\Figure2_RSkernel.png", replace


*_______________________________________________________________________________

*								MECHANISMS


* Preserve same sample used in the main analysis
xtset HHID y
replace sample_analysis=1 if F.sample_analysis==1 
keep if sample_analysis==1  


*-------------------------------------------------------------------------------

*** Allocation of per capita work time

	* Define control variables
		macro def CONT  c.hhsize c.age_head c.edu_adult c.dependency i.noadultmales i.incmlft ///
						c.ln_arableland c.cropdiv_SDI c.TLU i.intercrop ///
						i.wealth_index i.ssn_freemaize i.ssn_freefood i.ssn_FCfW i.ssn_childfeed 
		gen femaleperc = . 
	* Estimation
		eststo clear
		local activity "hhagr hhbus wage cas"
		foreach var of local activity {	
			replace femaleperc = Fp_`var'
			replace femaleperc = 0 if femaleperc==.
			xtreg PC_`var' ( $INDEP c.femaleperc )##$MM, fe vce(cluster EA) 
			estimate store PC_`var'
			xtreg PC_`var' ( $INDEP c.femaleperc $CONT )##$MM , fe vce(cluster EA) 
			estimate store PC_`var'_c
		}
	* Export Table 4
		estout PC_* , 	cells(b(star fmt(3)) se(par fmt(3))) starlevels( * 0.10 ** 0.05 *** 0.010) ///
			stat(N N_g r2 sigma_u sigma_e rho, fmt(0 0 2 2 2 2))  ///
			keep(1.drought 1.landmanagement1 1.drought#1.landmanagement1 femaleperc 1.drought#1.comm_marriagetype_matri 1.landmanagement1#1.comm_marriagetype_matri 1.drought#1.landmanagement1#1.comm_marriagetype_matri 1.comm_marriagetype_matri#c.femaleperc )
		esttab PC_* using  "$TAB/Table4_Mechanisms_worktime.tex",  /// 
			b(%9.3f) se(%9.3f) star(* 0.10 ** 0.05 *** 0.01) ///
			stat(N N_g r2 sigma_u sigma_e rho, fmt(0 0 2 2 2 2) label("N obs." "N. groups" "R^2"  "Sigma u" "Sigma e" "rho")) ///
			tex replace noconstant  nogaps nobaselevels noomitted ///
			coeflabels(1.drought "Drought" 1.landmanagement1 "Female Land Management" ///
					   1.drought#1.landmanagement1 "Drought x Female LM" 1.comm_marriagetype_matri "Norm (1=Matrilineal-Matrilocal)" ///
					    1.drought#1.comm_marriagetype_matri "Drought x Norm" /// 	
					   1.landmanagement1#1.comm_marriagetype_matri "Female LM x Norm" ///
					   1.drought#1.landmanagement1#1.comm_marriagetype_matri "Drought x Female LM x Norm" ///
					   femaleperc "Women's work (perc.)" 1.comm_marriagetype_matri#c.femaleperc "Women's work (perc.) x Norm" )  ///
			order(1.drought 1.landmanagement1 1.drought#1.landmanagement1 femaleperc 1.drought#1.comm_marriagetype_matri 1.landmanagement1#1.comm_marriagetype_matri 1.drought#1.landmanagement1#1.comm_marriagetype_matri 1.comm_marriagetype_matri#c.femaleperc ) ///
			keep(1.drought 1.landmanagement1 1.drought#1.landmanagement1 femaleperc 1.drought#1.comm_marriagetype_matri 1.landmanagement1#1.comm_marriagetype_matri 1.drought#1.landmanagement1#1.comm_marriagetype_matri 1.comm_marriagetype_matri#c.femaleperc) ///
			title("Table 4: Mechanism: allocation of per capita work time") 

*-------------------------------------------------------------------------------

*** Livestock holdings

	* Define sets of control variables
		macro def CONT  c.hhsize c.age_head c.edu_adult c.dependency  i.noadultmales i.incmlft ///
						c.ln_arableland c.cropdiv_SDI i.intercrop ///
						i.wealth_index i.ssn_freemaize i.ssn_freefood i.ssn_FCfW i.ssn_childfeed

		macro def CONTLV c.TLUsold c.TLUslaugh c.TLUbought c.TLUlost c.TLUgiveaway
			
	* Estimation
		eststo clear
		foreach var of varlist  TLU lvstkvalue_log { // 
				
				xtreg `var' ($INDEP $CONT )##$MM if year>2010, fe vce(cluster EA) 
				estimate store `var'
				xtreg `var' ($INDEP $CONT $CONTLV )##$MM if year>2010, fe vce(cluster EA) 
				estimate store `var'_c
				preserve
					replace `var'=. if TLU==0 
					xtreg `var' ($INDEP $CONT $CONTLV )##$MM if year>2010, fe vce(cluster EA) 
					estimate store `var'_pos
				restore
		}
	* Export Table 5
		estout  TLU* lvstkvalue_log*, cells(b(star fmt(3)) se(par fmt(3))) starlevels( * 0.10 ** 0.05 *** 0.010) ///
			stat(N N_g r2 sigma_u sigma_e rho, fmt(0 0 2 2 2 2)) ///
			keep(1.drought 1.landmanagement1 1.drought#1.landmanagement1 TLUsold TLUslaugh TLUbought TLUlost TLUgiveaway 1.drought#1.comm_marriagetype_matri 1.landmanagement1#1.comm_marriagetype_matri 1.drought#1.landmanagement1#1.comm_marriagetype_matri  1.comm_marriagetype_matri#c.TLUsold 1.comm_marriagetype_matri#c.TLUslaugh 1.comm_marriagetype_matri#c.TLUbought 1.comm_marriagetype_matri#c.TLUlost 1.comm_marriagetype_matri#c.TLUgiveaway)
		esttab 	TLU* lvstkvalue_log*	using  "$TAB/Table5_Mechanisms_livestock.tex",  /// 
				b(%9.3f) se(%9.3f) star(* 0.10 ** 0.05 *** 0.01) ///
				stat(N N_g r2 sigma_u sigma_e rho, fmt(0 0 2 2 2 2) ///
					 label("N obs." "N. groups" " $R^2$ "  " $\sigma_u$" " $\sigma_e$" " $\rho$")) ///
				tex replace noconstant  nogaps nobaselevels noomitted ///
				coeflabels( 1.drought "Drought" 1.landmanagement1 "Female Land Management" ///
						   1.drought#1.landmanagement1 "Drought x Female LM" 1.comm_marriagetype_matri "Norm (1=Matrilineal-Matrilocal)" ///
						   1.drought#1.comm_marriagetype_matri "Drought x Norm" /// 	
						   1.landmanagement1#1.comm_marriagetype_matri "Female LM x Norm" ///
						   1.drought#1.landmanagement1#1.comm_marriagetype_matri "Drought x Female LM x Norm" TLUsold "Livestock sold alive (TLU)" TLUslaugh "Livestock slaughtered (TLU)" TLUbought "Livestock bought (TLU)" TLUlost "Livestock lost or stolen (TLU)" TLUgiveaway "Livestock given away for free (TLU)" 1.comm_marriagetype_matri#c.TLUsold "Livestock sold alive (TLU) x Norm" 1.comm_marriagetype_matri#c.TLUslaugh  "Livestock slaughtered (TLU) x Norm" 1.comm_marriagetype_matri#c.TLUbought  "Livestock bought (TLU) x Norm" 1.comm_marriagetype_matri#c.TLUlost "Livestock lost or stolen (TLU) x Norm" 1.comm_marriagetype_matri#c.TLUgiveaway "Livestock given away (TLU) x Norm" )  ///
				order(1.drought 1.landmanagement1 1.drought#1.landmanagement1 TLUsold TLUslaugh TLUbought TLUlost TLUgiveaway ///
				   1.drought#1.comm_marriagetype_matri 1.landmanagement1#1.comm_marriagetype_matri 1.drought#1.landmanagement1#1.comm_marriagetype_matri  1.comm_marriagetype_matri#c.TLUsold 1.comm_marriagetype_matri#c.TLUslaugh 1.comm_marriagetype_matri#c.TLUbought 1.comm_marriagetype_matri#c.TLUlost 1.comm_marriagetype_matri#c.TLUgiveaway) ///
				keep(1.drought 1.landmanagement1 1.drought#1.landmanagement1 TLUsold TLUslaugh TLUbought TLUlost TLUgiveaway ///
				   1.drought#1.comm_marriagetype_matri 1.landmanagement1#1.comm_marriagetype_matri 1.drought#1.landmanagement1#1.comm_marriagetype_matri  1.comm_marriagetype_matri#c.TLUsold 1.comm_marriagetype_matri#c.TLUslaugh 1.comm_marriagetype_matri#c.TLUbought 1.comm_marriagetype_matri#c.TLUlost 1.comm_marriagetype_matri#c.TLUgiveaway) ///
				title("Table 5: Mechanism: livestock holdings") 



