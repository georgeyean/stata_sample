/*
 * Chow test
 */
use datasets/cfps_adult, clear
gen log_income=log(p_income+1)
drop if qq1101<0
drop if te4<0

// geneate var
gen male=cfps_gender
gen female=1-cfps_gender
tab te4, gen(edu)
local nulls "(male=female)"
local explans "male female"

foreach v of varlist qq1101 edu*{

    gen `v'f=`v'*female
    gen `v'm=`v'*male
    local nulls "`nulls' (`v'f=`v'm)"
    local explans "`explans' `v'f `v'm"
}
// regress 
reg log_income qq1101 edu1-edu7 if male==1
reg log_income qq1101 edu1-edu7 if female==1

// merge reg and check
reg log_income `explans', noconstant
di "`nulls'"
test `nulls'


/*
 * Model selection 
 */
clear
set obs 50
gen x=runiform()*3
gen y=exp(x)+rnormal()*2
local control ""
scalar K=0

forvalues i=1/10{

    gen x`i'=x^`i'
    local control "`control' x`i'"
    scalar K=K+1
    quietly: reg y `control'
    scalar r2=e(r2)
    scalar r2a=e(r2_a)
    scalar aic=-2*e(ll)+2*K
    scalar bic=-2*e(ll)+log(e(N))*K
    display `i' _skip r2 _skip r2a _skip aic _skip bic
}


/*
 * Wald OLS
 */
set seed 19880505
cap program drop dgp

program define dgp, rclass

	syntax [, obs(integer 1000) b(real -1.0)]
	drop _all
	set obs `obs'
	tempvar x y u
	tempname vkk V test_stat bk p
	gen `x'=rnormal()
	gen `u'=rnormal()^3
	gen `y'=`b'*`x'+`u'
	reg `y' `x' , robust
	local `bk'=_b[`x']
	mat `V'=e(V)
	local `vkk'=`V'[rownumb(`V',"`x'"),colnumb(`V',"`x'")]
	local `test_stat'=((``bk'')^2-1)^2/(4*(``bk'')^2*``vkk'')
	local `p'=1-chi2(1,``test_stat'')
	
	if ``p''<0.05{
		return scalar rejected=1
	}
	else{
		return scalar rejected=0
	}
	return scalar teststat=``test_stat''
end

// simulate
simulate rejected=r(rejected) teststat=r(teststat),reps(2000):dgp
su
hist teststat
simulate rejected=r(rejected) teststat=r(teststat),reps(2000):dgp, b(0)
su


/*
 * White test
 */
set seed 19880505
cap program drop dgp

program define dgp, rclass

	syntax [, obs(integer 100) robust b(real 0)]
	drop _all
	set obs `obs'
	tempvar x y sigma
	gen `x'=rnormal()
	gen `sigma'=sqrt(`x'^2)
	gen `y'=`b'*`x'+`sigma'*rnormal()
	quietly: reg `y' `x', `robust'
	return scalar b=_b[`x']
	return scalar se=_se[`x']
	
	if abs(_b[`x']/_se[`x'])>=1.96 {
		return scalar rejected=1
	}
	else {
		return scalar rejected=0
	}
end

simulate rejected=r(rejected) b=r(b) se=r(se),reps(1000):dgp
su
simulate rejected=r(rejected) b=r(b) se=r(se),reps(1000):dgp, robust
su


/*
 * multicolinearty test
 */
set seed 19880505
cap program drop dgp

program define dgp, rclass

	syntax [, obs(integer 20) b(real 0) mc(real 1)]
	drop _all
	set obs `obs'
	tempvar x1 x2 y sigma
	gen `x1'=rnormal()
	
	if `mc'==1{
		gen `x2'=3/sqrt(10)*`x1'+1/sqrt(10)*rnormal()
	}
	else {
		gen `x2'=rnormal()
	}
	
	gen `y'=`b'*`x1'+`x2'+rnormal()
	quietly: reg `y' `x1' `x2'
	return scalar b=_b[`x1']
	return scalar se=_se[`x1']
	
	if abs(_b[`x1']/_se[`x1'])>=invt(`obs'-3,0.975) {
		return scalar rejected=1
	}
	else {
		return scalar rejected=0
	}
	
	return scalar neg=_b[`x1']<0
end
// multicolinearty, size
simulate rejected=r(rejected) b=r(b) se=r(se) neg=r(neg),reps(2000):dgp
su
// multicolinearty, power
simulate rejected=r(rejected) b=r(b) se=r(se) neg=r(neg),reps(2000):dgp, b(1)
su
// no multicolinearty, power
simulate rejected=r(rejected) b=r(b) se=r(se) neg=r(neg),reps(2000):dgp, b(1) mc(0)
su
// multicolinearty, power with large N
simulate rejected=r(rejected) b=r(b) se=r(se) neg=r(neg),reps(2000):dgp, b(1) obs(100)
su
