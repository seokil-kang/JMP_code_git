estimation codes for "quantifying the fiscal backing for monetary policy"
Seokil Kang, 2021 sk86@iu.edu
-----------------------------------------------------------------------------------------------------
data_study.m: report data stats. It shares the exact code with data_importing function in data folder. It plots Figure 1, 12, 13
-----------------------------------------------------------------------------------------------------

-----------------------------------------------------------------------------------------------------
data: contains raw data
data_comparison.m: compare dataset with that of Cochrane(2021)
data_importing: import data for estimation for both DSGE & VAR
-----------------------------------------------------------------------------------------------------

-----------------------------------------------------------------------------------------------------
DSGE: contains DSGE estimation part
counterfactuals.m: Figure 4, 5
MP_shock_extracting.m: Figure 30
posterior_predictive.m: Figure 2, 15, 22 and Table 3, 5
posterior_predictive_fiscal.m: Figure 3
prior_predictive.m: Table 3
SMC_estimation_report.m: Table 2, Figure 14, 23 ~ 29

benchmark: contains model specifics
main_SMC.M: run the SMC estimation for the DSGE model
model: contains the structure of the model, particularly the file equilibrium_condition
functions: contains bunch of relevant submodules
prior: contains submodules prior density computation
gensys, sims_minwel: imported from Chris Sims

-----------------------------------------------------------------------------------------------------

-----------------------------------------------------------------------------------------------------
VAR: contains BVAR estimation part
proxy_iv_HD.m: Figure 6, 11, 17
proxy_iv_irf.m: Figure 7, 8, 9, 15, 16, 19, 20 Table 4, 6
sign_id_irf.m: Estimate with sign restriction.
proxy_iv_VD.m: Figure 10
BPSVAR: Bayesian Proxy SVAR(by Caldara&Herbst 2019)
main_BPSVAR.m: run the estimation and save results
IRF_DVD_BPSVAR.m: Figure 18, Table 7
HBVAR: VAR estimation for proxy-IV and sign restriction with Hierarchical BVAR(Giannone et al 2015)
mainfile_HBVAR.m: run the estimation and save results
SVAR: contains Bayesian SVAR estimation(Arias et al 2018)
mainfile.m: run the estimation and save results
report_result.m: Figure 21, Table 8