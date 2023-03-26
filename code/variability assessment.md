On page 16 in our manuscript, to assess the variability of quantiles, we consider two seed numbers.

AngCPD_Out=Ang_CPD(TPDM_P_hat = TPDM_Out$TPDM_P_hat, Add_c = 5, m = 85)
#[1] 0.2941503 1.4016199, 510 point masses
set.seed(237)
AngCPD_Out=Ang_CPD(TPDM_P_hat = TPDM_Out$TPDM_P_hat, Add_c = 6, m = 80)
#[1] [1] 0.3310329 1.4060339, 560 point masses

#stored data [1] 0.2780309 1.4152875, 459 pts
#second set [1] 0.3157301 1.4341370, 510 pts
#third set [1] 0.2811725 1.4058627, 560 pts
