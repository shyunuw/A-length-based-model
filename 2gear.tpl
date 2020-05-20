DATA_SECTION
  ////////////////////////////////////////////////////////
  // The number
  ////////////////////////////////////////////////////////
  init_int nages; // The number of age classes;
  init_int nlengths; // The number of length classes;
  init_int nmos; // The number of bi-months
  int r;
    !!r=1; //recruitment is defined as the pop size at age 1; 
  int ncohorts; // the number of chorts class;
    !!ncohorts=nmos; 
  
  ////////////////////////////////////////////////////////
  // Yield & CPUE data (Jigger & Purse-seine)
  ////////////////////////////////////////////////////////
  init_matrix YC(1,nmos,1,6); //fishery yield-CPUE ;
   vector YEAR(1,nmos); 
      !!YEAR=column(YC,1); 
    vector MONTH(1,nmos); 
      !!MONTH=column(YC,2); 
    vector JIG_CPUE(1,nmos); // Jigger;
      !!JIG_CPUE=column(YC,3);
    vector JIG_yield(1,nmos); 
      !!JIG_yield=column(YC,4);
    vector PS_CPUE(1,nmos); // Purse-seine 
      !!PS_CPUE=column(YC,5);
    vector PS_yield(1,nmos);
      !!PS_yield=column(YC,6);
  
  ////////////////////////////////////////////////////////
  // Derived Effort data.
  ////////////////////////////////////////////////////////
  vector JIG_effort(1,nmos);  //effort data derived
    !!JIG_effort=elem_div(JIG_yield,JIG_CPUE); 
  vector PS_effort(1,nmos);  
    !!PS_effort=elem_div(PS_yield,PS_CPUE);   
    
  ////////////////////////////////////////////////////////
  // Index for length frenquency data.
  ////////////////////////////////////////////////////////
  init_int indexMinmoLD; // index of the minimum bi-month for the length data; //2016. 5&6;
  int nmosLD;       //bi-months for the length data ;
    !!nmosLD=nmos-indexMinmoLD+1;
     
  ////////////////////////////////////////////////////////
  // Discrete length classes .
  ////////////////////////////////////////////////////////
  init_vector x(1,nlengths); //discrete lengths 0.5,11.5,...,33.5(cm); the # of x = 34 
  vector L(1,nlengths);      //discrete lengths 0.5,11.5,...,33.5(cm); the # of x = 34
    !!L=x;
  
  ////////////////////////////////////////////////////////
  // Length frequency data.
  ////////////////////////////////////////////////////////
  init_matrix JIG_lengthfrq(1,nmosLD,1,nlengths); // size matrix;
  init_matrix PS_lengthfrq(1,nmosLD,1,nlengths); // size matrix;
  
  ////////////////////////////////////////////////////////
  // Sample size for length data.
  ////////////////////////////////////////////////////////
  ivector JIG_SamSize(indexMinmoLD,nmos);   //bimonthly sample size for the JIG length data;
  ivector PS_SamSize(indexMinmoLD,nmos);   //bimonthly sample size for the PS length data; 
  
  ////////////////////////////////////////////////////////
  // delta & eta: For missing data cell in PS data.
  ////////////////////////////////////////////////////////
  vector PS_delta(1,nmos);
  vector PS_eta(1,nmos);
  
  ////////////////////////////////////////////////////////
  // Weight parameter for the each likelihood function
  ////////////////////////////////////////////////////////
  init_number lambda1; // weight in the objective function for composition_JIG;
  init_number lambda2; // weight in the objective function for composition_PS;
  init_number JIG_sig2_yield; // weight in the objective function for abundance_JIG (Yield);
  init_number PS_sig2_yield; // weight in the objective function for abundance_PS (Yield);
  init_number lambda3; // weight in the objective function for length-weight data;
  init_number lambda4; // weight in the objective function for Fecundity data;
    
  ////////////////////////////////////////////////////////
  // !!SWITCH_SECTION  !!
  ////////////////////////////////////////////////////////
     // 1. Switch [Stock-Recruitment model]
         init_int SR_switch;
         int nRecruits;
         int nProjection;
         
         !!nRecruits=indexMinmoLD;
         !!nProjection=1;
         !!if(SR_switch == 0) {
               !!nRecruits=nmos;
               !!nProjection=0;
         !!}
              // 1.Choice [Stock-Recruitment model];
                  //1-1. 
                  init_int choiceSR;
                  //1-2.
                   init_int Maturation_switch;
                   init_vector par_logistic(1,2); //For the maturation rate by length class; (Jo et al. 2019).;
                   vector maturation(1,nlengths); 
                   number b0;
                   number b1;
                   !!b0=0.0; 
                   !!b1=0.0; 
                   !!if(SR_switch == 1) {
                       !!if(Maturation_switch == 1) {
                            !!b0=par_logistic(1); 
                            !!b1=par_logistic(2); 
                            !!for(int xind=1;xind<=nlengths;xind++)
                                  !!maturation(xind)=mfexp(b0+b1*x(xind))/(1+mfexp(b0+b1*x(xind)));
                       !!}
                       !!else if(Maturation_switch == 0)
                           !!for(int xind=1;xind<=nlengths;xind++)
                                   !!maturation(xind)=1.0;
                   !!};
                  //1-2.
                  init_int Fecundity_switch;
                  init_int nFC;
                  init_vector lengthFC(1,nFC);
                  init_vector EggFC(1,nFC);

     // 2. Switch [Effective sample size]
         init_int EFF_Ssize_switch;
         init_int EFF_SamSize;
         
     // 3. Switch [Instantaneous natural mortality: M]
         init_int M_switch;
         init_number M_input; // M: (0.1/month) From T.kaga.
         
     // 4. Switch [Allometric relationship: AL]
         init_int LW_switch;
         init_int nAL; //the sample size for the AL
         init_matrix LW(1,2,1,nAL);  // length-weight data for AL
         vector lengthAL(1,nAL); // in cm
         vector weightAL(1,nAL); // in gram 
         !!lengthAL=LW(1);  //1st row
         !!weightAL=LW(2);  //2nd row
     // 5. Switch [Age Plus Group]
         init_int Plus_switch;
     // 6. Switch [Selection Phase number] // Note that if you declare phase number as '-1' for any parameter, ADMB doesn't estimate that parameter.
         int nPARAMETERS;
         init_int phase_mu_r;
         !!if(phase_mu_r != -1)
              !!nPARAMETERS+=1;
         init_int phase_log_sig2_r;
         !!if(SR_switch == 1 && choiceSR == 0 && Fecundity_switch ==1 )
              !!phase_log_sig2_r=1;
         !!if(phase_log_sig2_r != -1)
              !!nPARAMETERS+=1;
         init_int phase_Linf;
         !!if(phase_Linf != -1)
              !!nPARAMETERS+=1;
         init_int phase_log_kappa;
         !!if(phase_log_kappa != -1)
              !!nPARAMETERS+=1;
         init_int phase_sigmaL;
         !!if(phase_sigmaL != -1)
              !!nPARAMETERS+=1;
         init_int phase_JIG_log_q;
         !!if(phase_JIG_log_q != -1)
              !!nPARAMETERS+=1;
         init_int phase_JIG_log_gamma;
         !!if(phase_JIG_log_gamma != -1)
              !!nPARAMETERS+=1;
         init_int phase_JIG_log_L50;
         !!if(phase_JIG_log_L50 != -1)
              !!nPARAMETERS+=1;
         init_int phase_PS_log_q;
         !!if(phase_PS_log_q != -1)
              !!nPARAMETERS+=1;
         init_int phase_PS_log_gamma;
         !!if(phase_PS_log_gamma != -1)
              !!nPARAMETERS+=1;
         init_int phase_PS_log_L50;
         !!if(phase_PS_log_L50 != -1)
              !!nPARAMETERS+=1;
         init_int phase_logRec;
         !!if(phase_logRec != -1)
              !!nPARAMETERS+=1*nRecruits;
         init_int phase_log_Egg_M;
         !!if(SR_switch == 1 && choiceSR == 0 && phase_log_Egg_M != -1)
              !!nPARAMETERS+=1;
         init_int phase_log_aFC;
         init_int phase_log_bFC;
         !!if(SR_switch == 1 && choiceSR == 0 && Fecundity_switch == 0) {
             !!phase_log_aFC=-1;
             !!phase_log_bFC=-1;
         !!};
         !!if(SR_switch == 1 && choiceSR == 0 && phase_log_aFC != -1)
              !!nPARAMETERS+=1;
         !!if(SR_switch == 1 && choiceSR == 0 && phase_log_bFC != -1)
              !!nPARAMETERS+=1;
         init_int phase_logaSR;
         !!if(SR_switch == 1 && choiceSR != 0 && phase_logaSR != -1)
              !!nPARAMETERS+=1;
         init_int phase_logbSR;
         !!if(SR_switch == 1 && choiceSR != 0 && phase_logbSR != -1)
              !!nPARAMETERS+=1;
         init_int phase_Sp_rate;
         !!if(SR_switch == 1 && phase_Sp_rate != -1)
              !!nPARAMETERS+=1;
         init_int phase_log_parM;
         !!if(M_switch == 1 && phase_log_parM != -1)
              !!nPARAMETERS+=1;
         init_int phase_log_powerM;
         !!if(M_switch == 1 && phase_log_powerM != -1)
              !!nPARAMETERS+=1;
         init_int phase_log_aWL;
         init_int phase_bWL;
         !!if(LW_switch == 0) {
             !!phase_log_aWL=-1;
             !!phase_bWL=-1;
         !!};
         !!if(LW_switch == 1 && phase_log_aWL != -1)
              !!nPARAMETERS+=1;
         !!if(LW_switch == 1 && phase_bWL != -1)
              !!nPARAMETERS+=1;

  ////////////////////////////////////////////////////////
  // !! Choice a pin file !!
  ////////////////////////////////////////////////////////
         !!if(SR_switch == 0) {
               !!if(M_switch == 0) {
                    !!if(LW_switch == 0) 
                          !!ad_comm::change_pinfile_name("2gear.pin");
                    !!else if(LW_switch == 1) 
                          !!ad_comm::change_pinfile_name("2gear_LW.pin");
               !!}
               !!if(M_switch == 1) {
                    !!if(LW_switch == 0) 
                          !!ad_comm::change_pinfile_name("2gear_M.pin");
                    !!else if(LW_switch == 1) 
                          !!ad_comm::change_pinfile_name("2gear_M_LW.pin");
               !!}
         !!}
         !!else if(SR_switch == 1) {
               !!if(choiceSR == 0) {
                    !!if(M_switch == 0) {
                          !!if(LW_switch == 0) 
                                !!ad_comm::change_pinfile_name("Egg.pin");
                          !!else if(LW_switch == 1) 
                                !!ad_comm::change_pinfile_name("Egg_LW.pin");
                    !!}
                    !!else if(M_switch == 1)
                          !!if(LW_switch == 0) 
                                !!ad_comm::change_pinfile_name("Egg_M.pin");
                          !!else if(LW_switch == 1) 
                                !!ad_comm::change_pinfile_name("Egg_M_LW.pin");
               !!}
               !!else if(choiceSR == 1) {
                    !!if(M_switch == 0) {
                          !!if(LW_switch == 0) 
                                !!ad_comm::change_pinfile_name("B_H.pin");
                          !!else if(LW_switch == 1) 
                                !!ad_comm::change_pinfile_name("B_H_LW.pin");
                    !!}
                    !!else if(M_switch == 1)
                          !!if(LW_switch == 0) 
                                !!ad_comm::change_pinfile_name("B_H_M.pin");
                          !!else if(LW_switch == 1) 
                                !!ad_comm::change_pinfile_name("B_H_M_LW.pin");
               !!}
               !!else if(choiceSR == 2) {
                    !!if(M_switch == 0) {
                          !!if(LW_switch == 0) 
                                !!ad_comm::change_pinfile_name("Ricker.pin");
                          !!else if(LW_switch == 1) 
                                !!ad_comm::change_pinfile_name("Ricker_LW.pin");
                    !!};
                    !!if(M_switch == 1)
                          !!if(LW_switch == 0) 
                                !!ad_comm::change_pinfile_name("Ricker_M.pin");
                          !!else if(LW_switch == 1) 
                                !!ad_comm::change_pinfile_name("Ricker_M_LW.pin");
               !!};
         !!};
  ////////////////////////////////////////////////////////
  // !! Check the Input values !!
  ////////////////////////////////////////////////////////
  !!cout<<"nages: "<<nages<<endl;
  !!cout<<"nlengths: "<<nlengths<<endl;
  !!cout<<"nmos: "<<nmos<<endl;
  !!cout<<"JIG_CPUE: "<<JIG_CPUE<<endl;
  !!cout<<"JIG_yield: "<<JIG_yield<<endl;
  !!cout<<"PS_CPUE: "<<PS_CPUE<<endl;
  !!cout<<"PS_yield: "<<PS_yield<<endl;
  !!cout<<"JIG_effort: "<<JIG_effort<<endl;
  !!cout<<"PS_effort: "<<PS_effort<<endl;
  !!cout<<"indexMinmoLD: "<<indexMinmoLD<<endl;
  !!cout<<"nmosLD: "<<nmosLD<<endl;
  !!cout<<"x: "<<x<<endl;
  !!cout<<"L: "<<L<<endl;
  !!cout<<"JIG_lengthfrq: "<<JIG_lengthfrq<<endl;
  !!cout<<"PS_lengthfrq: "<<PS_lengthfrq<<endl;
  !!cout<<"lambda1: "<<lambda1<<endl;
  !!cout<<"lambda2: "<<lambda2<<endl;
  !!cout<<"JIG_sig2_yield: "<<JIG_sig2_yield<<endl;
  !!cout<<"PS_sig2_yield: "<<PS_sig2_yield<<endl;
  !!cout<<"lambda3: "<<lambda3<<endl;
  !!cout<<"lambda4: "<<lambda4<<endl;
  !!cout<<"SR_switch: "<<SR_switch<<endl;
  !!cout<<"choiceSR: "<<choiceSR<<endl;
  !!cout<<"nFC: "<<nFC<<endl;
  !!cout<<"lengthFC: "<<lengthFC<<endl;
  !!cout<<"EggFC: "<<EggFC<<endl;  
  !!cout<<"Maturation_switch: "<<Maturation_switch<<endl;
  !!cout<<"par_logistic: "<<par_logistic<<endl;
  !!cout<<"EFF_Ssize_switch: "<<EFF_Ssize_switch<<endl;
  !!cout<<"EFF_SamSize: "<<EFF_SamSize<<endl;
  !!cout<<"M_switch: "<<M_switch<<endl;
  !!cout<<"M_input: "<<M_input<<endl;
  !!cout<<"LW_switch: "<<LW_switch<<endl;
  !!cout<<"nAL: "<<nAL<<endl;
  !!cout<<"LW: "<<LW<<endl;
  !!cout<<"lengthAL: "<<lengthAL<<endl;
  !!cout<<"weightAL: "<<weightAL<<endl;
  !!cout<<"Plus_switch: "<<Plus_switch<<endl;
  !!cout<<"phase_mu_r: "<<phase_mu_r<<endl;
  !!cout<<"phase_log_sig2_r: "<<phase_log_sig2_r<<endl;
  !!cout<<"phase_Linf: "<<phase_Linf<<endl;
  !!cout<<"phase_log_kappa: "<<phase_log_kappa<<endl;
  !!cout<<"phase_sigmaL: "<<phase_sigmaL<<endl;
  !!cout<<"phase_JIG_log_q: "<<phase_JIG_log_q<<endl;
  !!cout<<"phase_JIG_log_gamma: "<<phase_JIG_log_gamma<<endl;
  !!cout<<"phase_JIG_log_L50: "<<phase_JIG_log_L50<<endl;
  !!cout<<"phase_PS_log_q: "<<phase_PS_log_q<<endl;
  !!cout<<"phase_PS_log_gamma: "<<phase_PS_log_gamma<<endl;
  !!cout<<"phase_logRec: "<<phase_logRec<<endl;
  !!cout<<"phase_log_Egg_M: "<<phase_log_Egg_M<<endl;
  !!cout<<"phase_log_aFC: "<<phase_log_aFC<<endl;
  !!cout<<"phase_log_bFC: "<<phase_log_bFC<<endl;
  !!cout<<"phase_logaSR: "<<phase_logaSR<<endl;
  !!cout<<"phase_logbSR: "<<phase_logbSR<<endl;
  !!cout<<"phase_Sp_rate: "<<phase_Sp_rate<<endl;
  !!cout<<"phase_log_parM: "<<phase_log_parM<<endl;
  !!cout<<"phase_log_powerM: "<<phase_log_powerM<<endl;
  !!cout<<"phase_log_aWL: "<<phase_log_aWL<<endl;
  !!cout<<"phase_bWL: "<<phase_bWL<<endl;
  //!!exit(34);
PARAMETER_SECTION
  ////////////////////////////////////////////////////////
  // Free parameters
  ////////////////////////////////////////////////////////
  init_bounded_number mu_r(1.0,2.0,phase_mu_r);
  init_bounded_number log_sig2_r(-3.0,-1.2,phase_log_sig2_r);
  init_bounded_number Linf(33.0,40.0,phase_Linf);
  init_bounded_number log_kappa(-1.0,0.0,phase_log_kappa);   //growth para.;
  init_bounded_number sigmaL(0.0,0.1,phase_sigmaL);     //the uncertainty in the L_{a+1} equation;
  init_bounded_number JIG_log_q(-20.0,-12.0,phase_JIG_log_q);
  init_number JIG_log_gamma(phase_JIG_log_gamma);
  init_bounded_number JIG_log_L50(2.3,3.3,phase_JIG_log_L50);
  init_bounded_number PS_log_q(-20.0,-10.0,phase_PS_log_q);
  init_number PS_log_gamma(phase_PS_log_gamma);
  init_bounded_number PS_log_L50(2.3,3.3,phase_PS_log_L50);
  init_bounded_vector logRec(1,nRecruits,0.0,35.0,phase_logRec);  // log Recruits(number=6);
  !!cout<<"mu_r: "<<mu_r<<endl;
  !!cout<<"log_sig2_r: "<<log_sig2_r<<endl;
  !!cout<<"Linf: "<<Linf<<endl;
  !!cout<<"log_kappa: "<<log_kappa<<endl;
  !!cout<<"sigmaL: "<<sigmaL<<endl;
  !!cout<<"JIG_log_q: "<<JIG_log_q<<endl;
  !!cout<<"JIG_log_gamma: "<<JIG_log_gamma<<endl;
  !!cout<<"JIG_log_L50: "<<JIG_log_L50<<endl;
  !!cout<<"PS_log_q: "<<PS_log_q<<endl;
  !!cout<<"PS_log_gamma: "<<PS_log_gamma<<endl;
  !!cout<<"PS_log_L50: "<<PS_log_L50<<endl;
  !!cout<<"logRec: "<<logRec<<endl;
  
  !!if(SR_switch == 1) {
       !!if(choiceSR == 0) {
           init_number log_Egg_M(phase_log_Egg_M);
           init_number log_aFC(phase_log_aFC);
           init_number log_bFC(phase_log_bFC);
           number Egg_M;
           vector Eggs(indexMinmoLD,nmos+nProjection);
           vector fecundity(1,nlengths);
           number aFC;
           number bFC;
           number sig2FC; 
           vector E_FC(1,nFC);
           !!cout<<"log_Egg_M: "<<log_Egg_M<<endl;
           !!cout<<"log_aFC: "<<log_aFC<<endl;
           !!cout<<"log_bFC: "<<log_bFC<<endl;
       !!}
       !!else if(choiceSR != 0) {
           init_number logaSR(phase_logaSR);
           init_number logbSR(phase_logbSR);
           !!cout<<"logaSR: "<<logaSR<<endl;
           !!cout<<"logbSR: "<<logbSR<<endl;
       !!};
       init_bounded_number Sp_rate(0.0,1.0,phase_Sp_rate);
       !!cout<<"Sp_rate: "<<Sp_rate<<endl;
       number aSR;
       number bSR;
       3darray SpawnersL(1,nages,1,ncohorts+nProjection,1,nlengths);
       matrix Spawners(1,nages,1,nmos+nProjection);
  !!};
  !!if(M_switch == 1) {
       init_bounded_number log_parM(0.0,10.5,phase_log_parM); //2.5,10.5
       !!cout<<"log_parM: "<<log_parM<<endl;
       number parM;
       init_number log_powerM(phase_log_powerM);  //power (=beta*b1) in M = par*W^(-1.0*power)
       !!cout<<"log_powerM: "<<log_powerM<<endl;
       number powerM;
                    //ad hoc value = 0.305, where 0.305 from Lorenzen (1996);
  !!};
  init_number log_aWL(phase_log_aWL);
  number aWL; 
  init_bounded_number bWL(2.0,3.3,phase_bWL);
  !!cout<<"log_aWL: "<<log_aWL<<endl;
  !!cout<<"bWL: "<<bWL<<endl;
  vector Wt(1,nlengths);
  number sig2AL; 
  vector E_W(1,nAL);

  
  ////////////////////////////////////////////////////////
  // Derived quantites
  ////////////////////////////////////////////////////////
  // 1. Exponential form of free parameters
  vector Recruits(1,nmos+nProjection);
  number JIG_q; 
  number JIG_gamma; 
  number JIG_L50;
  number PS_q; 
  number PS_gamma;
  number PS_L50;
  number sig2_r;
  number kappa;
  
  // 2. Growth & Length distribution
  number Rho; // EXP[-kappa];
  number kkk; //for the cumulative purpose;
  matrix f(1,nages,1,nlengths);  //length frequency as pmf;
  3darray pp(2,nages,1,nlengths,1,nlengths); //pp(1 To Ages,1 To x,1 To L);
                                          // x --> (growth) --> L
  vector Mu(1,nlengths); //differ by length class;
  vector SS(1,nages); //assumed to be constant over length classes; indp. length;
  vector Mean_logN_L(2,nages);
  vector Var_logN_L(2,nages);
  
  // 3. Motality
  vector M(1,nlengths);
  vector JIG_Sel(1,nlengths);
  vector PS_Sel(1,nlengths);
  vector JIG_F_mo(1,nmos); 
  vector PS_F_mo(1,nmos);
  matrix JIG_F_tx(1,nmos,1,nlengths);
  matrix PS_F_tx(1,nmos,1,nlengths);
  matrix Z(1,nmos,1,nlengths); // Total Mortality;
  matrix ExpZ(1,nmos+nProjection,1,nlengths);  

  // 3-1. Catch
  number JIG_CNum;
  number PS_CNum;
  number JIG_CWt;
  number PS_CWt;
  vector JIG_Yieldhat(1,nmos);
  vector PS_Yieldhat(1,nmos);
  matrix JIG_Catch(1,nmos,1,nlengths);
  matrix PS_Catch(1,nmos,1,nlengths);
  3darray JIG_ENx(1,nages,1,ncohorts,1,nlengths); 
  vector JIG_EN(1,ncohorts); 
  vector JIG_EB(1,ncohorts); 
  3darray PS_ENx(1,nages,1,ncohorts,1,nlengths); 
  vector PS_EN(1,ncohorts); 
  vector PS_EB(1,ncohorts);
  
  // 4. Abundance 
  3darray NL(1,nages,1,ncohorts+nProjection,1,nlengths); 
  matrix N(1,nages,1,nmos+nProjection);  //dimension index is different from the VB code; 
  number SumP;
  vector p(1,nlengths); 
  vector Pop(1,nmos+nProjection); 
  vector B(1,ncohorts+nProjection);
  !!if(Plus_switch == 1) 
     vector N_plus(indexMinmoLD,nmos+nProjection);
  
  // 5. To Calculate Objective Function .
  number maxloglike; // for calculation of AIC
  
  // 5-1. part1: Catch_time(Length) ~ multinomial distribution
  matrix JIG_p_hat(indexMinmoLD,nmos,1,nlengths);
  matrix PS_p_hat(indexMinmoLD,nmos,1,nlengths);
  matrix JIG_LF(indexMinmoLD,nmos,1,nlengths);    
  matrix PS_LF(indexMinmoLD,nmos,1,nlengths);
  matrix EFF_JIG_lengthfrq(1,nmos-indexMinmoLD+1,1,nlengths);
  matrix EFF_PS_lengthfrq(1,nmos-indexMinmoLD+1,1,nlengths);
  number JIG_logmult; //log(multinomial);
  number PS_logmult;
  
  // 5-2. part2: log(yield) ~ Normal Distribution
  vector JIG_elem_obj2(indexMinmoLD,nmos);  //elements in part 2 of the objective function;
  vector PS_elem_obj2(indexMinmoLD,nmos);
  number lognormal;  //log(normal); 
  
  // 5-3. part3: only exist when turning on 'LW_switch'
  // Weight ~ Normal Distribution
  // Todarodes pacificus's L-W data doesn't appear goodness of fit when considering uncertainty as multiplicative error.
  // look at the 'Line 375 ~ 382'
  
  // 6. Objective function
  objective_function_value Q;  //negative loglikelihood
  
  // 6. AIC
  number aic;
  
PRELIMINARY_CALCS_SECTION
  //re-arrangement of data for calculation purposes;
  //bi-monthly sample size for the length data
  int i;
  JIG_SamSize=0;
  PS_SamSize=0;
  for(int m=indexMinmoLD;m<=nmos;m++) {
      i=m-indexMinmoLD+1;
      JIG_SamSize(m)=sum(JIG_lengthfrq(i));
      PS_SamSize(m)=sum(PS_lengthfrq(i));
  };
  // delta & eta for missing data cell.
  for(int i=1;i<=nmos;i++) { 
      if(PS_effort(i)==0)
          PS_delta(i)=PS_effort(i)+1.0;
      else
          PS_delta(i)=PS_effort(i)*0.0;
  };
  //cout<<"PS_delta : "<<PS_delta<<endl;
  for(int i=1;i<=nmos;i++) { 
      if(PS_delta(i)==1.0)
          PS_eta(i)=0.0;
      else
          PS_eta(i)=1.0;
  };
  //exit(34);
PROCEDURE_SECTION 
  JIG_gamma=mfexp(JIG_log_gamma);
  PS_gamma=mfexp(PS_log_gamma);
  JIG_L50=mfexp(JIG_log_L50);
  PS_L50=mfexp(PS_log_L50);
  JIG_q=mfexp(JIG_log_q);
  PS_q=mfexp(PS_log_q);
  sig2_r=mfexp(log_sig2_r);
  kappa=mfexp(log_kappa);
  if(SR_switch == 1) 
     if(choiceSR != 0) {
         aSR=mfexp(logaSR);
         bSR=mfexp(logbSR);
     };
  
  Q=0.0;
  
  aWL=mfexp(log_aWL); 
  Wt=aWL*pow(x,bWL)/1000;  //gram to(=>) kg
  
  //instantaneous natural mortality by length class
  if(M_switch == 0) 
       for(int xind=1; xind<=nlengths; xind++) 
            M(xind)=M_input; 
  else if(M_switch == 1) {
       parM=mfexp(log_parM);   
       powerM=mfexp(log_powerM);   
       for(int xind=1; xind<=nlengths; xind++) 
            M(xind)=parM*pow(aWL,bWL)*pow(x(xind),-1.0*powerM*bWL);
  };
  cout<<"debug_1: "<<endl;
  
  Recruits=0.0;
  for(int i=1;i<=nRecruits;i++)
      Recruits(i)=mfexp(logRec(i));
      
  cout<<"debug_2: "<<endl;
   
  //Calculate the gear selectivity for each length class (JIG & PS & DRAG); 
  JIG_Sel=1.0/(1.0+mfexp(-1.0*JIG_gamma*(x-JIG_L50)));
  PS_Sel=1.0/(1.0+mfexp(-1.0*PS_gamma*(x-PS_L50)));
  
  cout<<"debug_3: "<<endl;
  for(int t=1;t<=nmos;t++) 
     for(int xind=1;xind<=nlengths;xind++) { 
         JIG_F_mo(t)=JIG_q*JIG_effort(t);
         PS_F_mo(t)=PS_q*PS_effort(t);
         JIG_F_tx(t,xind)=JIG_F_mo(t)*JIG_Sel(xind); //JIG_F_tx; 
         PS_F_tx(t,xind)=PS_F_mo(t)*PS_Sel(xind); //PS_F_tx;
         Z(t,xind)=M(xind)+JIG_F_tx(t,xind)+PS_F_tx(t,xind);
         ExpZ(t,xind)=mfexp(-1.0*Z(t,xind));
     };
     
  cout<<"debug_4: "<<endl;
     
  //Calculating the length frequency dsn of the recruits;
  //recruitment is at 'TWO months' of age;
  //f(1,x) is the same for all cohorts so this can be in the initial calculations;
  Rho=mfexp(-1.0*kappa);
  SS(1)=sig2_r; 
  kkk=0.0;

  for(int xind=1;xind<=nlengths;xind++) {
     f(1,xind)=mfexp(-1.0*square(x(xind)-mu_r)/(2.0*SS(1)));
     kkk=kkk+f(1,xind);      //  page533, Xi_r (Quinn.);
  };
   
  for(int xind=1;xind<=nlengths;xind++) {
      f(1,xind)=f(1,xind)/kkk;  //normalize;  //
      Mu(xind)=Linf*(pow((x(xind)/Linf),Rho))*mfexp(pow(sigmaL,2.0)/2.0);
  };
  
  for(int a=2;a<=nages;a++) {
    Mean_logN_L(a)= (1.0-pow(Rho,a-1))*log(Linf)+pow(Rho,a-1)*(log(mu_r)-0.5*(sig2_r/pow(mu_r,2.0)));
    Var_logN_L(a)=pow(Rho,2.0*(a-1))*(sig2_r/pow(mu_r,2.0))+pow(sigmaL,2.0)*((1-pow(Rho,(2.0*(a-1))))/(1-pow(Rho,2.0)));
  };
  for(int a=2;a<=nages;a++)
      SS(a)=mfexp(2.0*(Mean_logN_L(a)+Var_logN_L(a)))-mfexp(2.0*Mean_logN_L(a)+Var_logN_L(a));
  
  for(int a=2;a<=nages;a++) 
     for(int xind=1;xind<=nlengths;xind++)  {
         kkk=0.0; // normalizing constant
         for(int Lind=1;Lind<=nlengths;Lind++)  {
             pp(a,Lind,xind)=0.0;
             if(Lind >= xind) {
                   pp(a,Lind,xind)=mfexp(-1.0/(2.0*SS(a))*square(L(Lind)-Mu(xind)));
                   
                  kkk=kkk+pp(a,Lind,xind);
             };
         };
         for(int Lind=1;Lind<=nlengths;Lind++)
             pp(a,Lind,xind)=pp(a,Lind,xind)/kkk;
     };
  
  cout<<"debug_5: "<<endl;
  //Start of cohort loop
  int m;
  N=0.0;
  Spawners=0.0;
  for(int c=1;c<=ncohorts+nProjection;c++)  {
     int a=1; 
     m=c;        //m: time index
     N(a,m)=0.0;
     if(SR_switch == 0) {
         for(int xind=1;xind<=nlengths;xind++)  {
                NL(a,m,xind)=0.0;
                NL(a,m,xind)=Recruits(m)*f(a,xind);   // recruit number (2015.07&08 - 2018.1&12)
                N(a,m)=N(a,m)+NL(a,m,xind);    //Be careful that N has ages x months; //different from the VB code
          };
     }
     else if(SR_switch == 1) {
         if(c <= nRecruits) {
            for(int xind=1;xind<=nlengths;xind++)  {
                NL(a,m,xind)=0.0;
                NL(a,m,xind)=Recruits(m)*f(a,xind);   // recruit number (2015.07&08 - 2018.1&12)
                N(a,m)=N(a,m)+NL(a,m,xind);    //Be careful that N has ages x months; //different from the VB code
            };
         }
         else {
               if(choiceSR == 0) {
                   Egg_M=mfexp(log_Egg_M);
                   aFC=mfexp(log_aFC);
                   bFC=mfexp(log_bFC);
                   fecundity=aFC*pow(x,bFC);
                   Eggs(m-1)=0.0;
                   for(int xind=1;xind<=nlengths;xind++) 
                       Eggs(m-1)+=(SpawnersL(5,m-1,xind)+SpawnersL(6,m-1,xind))*fecundity(xind);
                   N(a,m)=Eggs(m-1)*mfexp(-1.0*Egg_M)    ;  //Eggs
               };
               if(choiceSR == 1) 
                   N(a,m)=sum(column(Spawners,m-1))/(aSR+bSR*sum(column(Spawners,m-1)));  //B-H   
               if(choiceSR == 2) 
                   N(a,m)=aSR*sum(column(Spawners,m-1))*mfexp(-1.0*bSR*sum(column(Spawners,m-1)));  //Ricker   
                   
                for(int xind=1;xind<=nlengths;xind++) {
                      NL(a,m,xind)=0.0;
                      NL(1,c,xind) = N(1,c)*f(1,xind);
                 };
         };
     };
      for(int a=2;a<=nages;a++)   { //note that a starts at 2, not 1; 
         m=a+c-1;
         if(m<=nmos+nProjection) {
               SumP=0.0; 
               for(int Lind=1;Lind<=nlengths;Lind++)  {
                  p(Lind)=0;    
                  for(int xind=1;xind<=nlengths;xind++)
                        p(Lind)=p(Lind)+f(a-1,xind)*ExpZ(m-1,xind)*pp(a,Lind,xind); 
                         
                  SumP=SumP+p(Lind);
               };
               if(SR_switch == 0) {
                  N(a,m)=0.0;
                  for(int Lind=1;Lind<=nlengths;Lind++)  {
                        NL(a,m,Lind)=0.0;
                        f(a,Lind)=p(Lind)/SumP;  //normalize;
                        NL(a,m,Lind)=N(a-1,m-1)*p(Lind);
                        N(a,m)=N(a,m)+NL(a,m,Lind);
                  };
                  if(Plus_switch == 1 && a == nages) {
                       N_plus(m)=0.0;
                       for(int Lind=1;Lind<=nlengths;Lind++)  
                            for(int xind=1;xind<=nlengths;xind++)  {
                                  NL(a,m,Lind)+=NL(a,m-1,xind)*ExpZ(m-1,xind)*pp(a,Lind,xind);
                                  N(a,m)+=NL(a,m-1,xind)*ExpZ(m-1,xind)*pp(a,Lind,xind);
                                  N_plus(m)+=NL(a,m-1,xind)*ExpZ(m-1,xind)*pp(a,Lind,xind);
                            };
                  };
               }
               else if(SR_switch == 1) {
                  if(a!=nages) {
                        N(a,m)=0.0;
                        for(int Lind=1;Lind<=nlengths;Lind++)  {
                              NL(a,m,Lind)=0.0;
                              f(a,Lind)=p(Lind)/SumP;  //normalize;
                              NL(a,m,Lind)=N(a-1,m-1)*p(Lind);
                              N(a,m)=N(a,m)+NL(a,m,Lind);
                        };
                  }
                  else if(a == nages) {
                        Spawners(a-1,m-1)=0.0;
                        for(int Lind=1;Lind<=nlengths;Lind++)  {
                              SpawnersL(a-1,m-1,Lind)=0.0;
                              if(Maturation_switch == 0)
                                  SpawnersL(a-1,m-1,Lind)=NL(a-1,m-1,Lind)*Sp_rate;
                              else if(Maturation_switch == 1)
                                  SpawnersL(a-1,m-1,Lind)=NL(a-1,m-1,Lind)*Sp_rate*maturation(Lind);
                                  
                              Spawners(a-1,m-1)=Spawners(a-1,m-1)+SpawnersL(a-1,m-1,Lind); 
                        };
                        f(a)=0.0;
                        for(int xind=1;xind<=nlengths;xind++)  {
                              if(Maturation_switch==0)
                                    f(a,xind)=f(a-1,xind)*ExpZ(m-1,xind)*(1.0-Sp_rate); //Note that: f(a,xind) is age 5 at the end of time (m-1)
                              else if(Maturation_switch==1)
                                    f(a,xind)=f(a-1,xind)*ExpZ(m-1,xind)*(1.0-Sp_rate*maturation(xind)); //Note that: f(a,xind) is age 5 at the end of time (m-1)
                        };
                        for(int Lind=1;Lind<=nlengths;Lind++)  {
                              p(Lind)=0;
                              for(int xind=1;xind<=nlengths;xind++) {
                                    p(Lind)=p(Lind)+f(a,xind)*pp(a,Lind,xind);  //Note that: f(a,xind) is age 5 at the end of time (m-1)
                              };
                        };
                        N(a,m)=0.0;
                        Spawners(a,m)=0.0;
                        for(int Lind=1;Lind<=nlengths;Lind++)  {
                              NL(a,m,Lind)=0.0;
                              SpawnersL(a,m,Lind)=0.0;
                              NL(a,m,Lind)=N(a-1,m-1)*p(Lind);
                              if(Plus_switch == 1 && Maturation_switch == 1) {
                                    N_plus(m)=0.0;
                                    for(int xind=1;xind<=nlengths;xind++) {
                                                NL(a,m,Lind)+=NL(a,m-1,xind)*ExpZ(m-1,xind)*(1.0-maturation(xind))*pp(a,Lind,xind);
                                                N_plus(m)+=NL(a,m-1,xind)*ExpZ(m-1,xind)*(1.0-maturation(xind))*pp(a,Lind,xind);
                                    };
                              };
                              N(a,m)=N(a,m)+NL(a,m,Lind);
                              if(Maturation_switch == 0)
                                   SpawnersL(a,m,Lind)=NL(a,m,Lind); 
                              else if(Maturation_switch == 1)
                                   SpawnersL(a,m,Lind)=NL(a,m,Lind)*maturation(Lind); 
                                 
                              Spawners(a,m)=Spawners(a,m)+SpawnersL(a,m,Lind); 
                        };
                        f(a)=NL(a,m)/N(a,m);//Note that: f(a,xind) is age 6 at the begining of time (m)  //normalize; f(6) is a vector
                  };
               };
          };
      };
  };
  // Projection;
  if(SR_switch == 1) {
      for(int Lind=1;Lind<=nlengths;Lind++)  {
            SpawnersL(nages-1,nmos+nProjection,Lind)=0.0;
            if(Maturation_switch==0)
                 SpawnersL(nages-1,nmos+nProjection,Lind)=NL(nages-1,nmos+nProjection,Lind)*Sp_rate;
            else if(Maturation_switch==1)
                  SpawnersL(nages-1,nmos+nProjection,Lind)=NL(nages-1,nmos+nProjection,Lind)*Sp_rate*maturation(Lind);
            Spawners(nages-1,nmos+nProjection)=Spawners(nages-1,nmos+nProjection)+SpawnersL(nages-1,nmos+nProjection,Lind); 
      };
      if(choiceSR == 0){
           Eggs(nmos+nProjection)=0.0;
           for(int Lind=1;Lind<=nlengths;Lind++)  
                 Eggs(nmos+nProjection)+=(SpawnersL(nages-1,nmos+nProjection,Lind)+SpawnersL(nages,nmos+nProjection,Lind))*fecundity(Lind);
      };
  };
  cout<<"debug_6: "<<endl;
  
  if(SR_switch == 1)
      for(int i=indexMinmoLD+1;i<=nmos+nProjection;i++) 
            Recruits(i)=N(1,i);
            
  cout<<"debug_7: "<<endl;
  
  //note m starts at indexMinmoLD(=6) i.e, 2016yr May&June;
  //Sum through the population matrix to calculate;
  //expected catch distribution;
  //summing for each time by length class for all ages;
  for(int m=indexMinmoLD;m<=nmos;m++) { 
      JIG_CNum=0.0;
      JIG_CWt=0.0;
      JIG_Catch(m)=0.0;
      JIG_Yieldhat(m)=0.0;
      PS_CNum=0.0;
      PS_CWt=0.0;
      PS_Catch(m)=0.0;
      PS_Yieldhat(m)=0.0;
      Pop(m)=0.0;     
      B(m)=0.0;
      JIG_EN(m)=0.0;
      PS_EN(m)=0.0;
      JIG_EB(m)=0.0;
      PS_EB(m)=0.0;
      for(int a=1;a<=nages;a++) {
         for(int xind=1;xind<=nlengths;xind++) {
            JIG_CNum=(NL(a,m,xind))*(JIG_F_tx(m,xind)/Z(m,xind))*(1-ExpZ(m,xind));
            PS_CNum=(NL(a,m,xind))*(PS_F_tx(m,xind)/Z(m,xind))*(1-ExpZ(m,xind));
            JIG_ENx(a,m,xind)=NL(a,m,xind)*JIG_Sel(xind); //Exploitable population;
            PS_ENx(a,m,xind)=NL(a,m,xind)*PS_Sel(xind); //Exploitable population;
            JIG_EN(m)=JIG_EN(m)+JIG_ENx(a,m,xind);            //Exploitable population; 
            PS_EN(m)=PS_EN(m)+PS_ENx(a,m,xind);            //Exploitable population; 
            
            JIG_Catch(m,xind)=JIG_Catch(m,xind)+JIG_CNum;
            PS_Catch(m,xind)=PS_Catch(m,xind)+PS_CNum;
            JIG_CWt=JIG_CNum*Wt(xind);  //in kg
            PS_CWt=PS_CNum*Wt(xind);  
            B(m)=B(m)+NL(a,m,xind)*Wt(xind);
            JIG_EB(m)=JIG_EB(m)+JIG_ENx(a,m,xind)*Wt(xind);   //expoitable biomass; 
            PS_EB(m)=PS_EB(m)+PS_ENx(a,m,xind)*Wt(xind);   //expoitable biomass; 
            JIG_Yieldhat(m)=JIG_Yieldhat(m)+JIG_CWt; //in kg
            PS_Yieldhat(m)=PS_Yieldhat(m)+PS_CWt; //in kg
            
            Pop(m)=Pop(m)+NL(a,m,xind);   //population;
        };
     };
  };
  cout<<"debug_8: "<<endl;
  // Projection;
  if(SR_switch == 1) {
     Pop(nmos+nProjection)=0.0;
     B(nmos+nProjection)=0.0;
     for(int a=1;a<=nages;a++) 
          for(int xind=1;xind<=nlengths;xind++) {
                 Pop(nmos+nProjection)=Pop(nmos+nProjection)+NL(a,nmos+nProjection,xind);   //population;
                 B(nmos+nProjection)=B(nmos+nProjection)+NL(a,nmos+nProjection,xind)*Wt(xind);
          };
  }
  cout<<"debug_9: "<<endl;
  //exit(34);
  //part 1 of the objective funcion: multinomial;
  //The expected length-frequency 
  for(int m=indexMinmoLD;m<=nmos;m++)  { //note m starts at indexMinmoLD(=6);  
     for(int xind=1;xind<=nlengths;xind++)  {
         JIG_p_hat(m,xind)=JIG_Catch(m,xind)/sum(JIG_Catch(m));
         PS_p_hat(m,xind)=(PS_Catch(m,xind)+PS_delta(m))/(sum(PS_Catch(m)+PS_delta(m)));
         JIG_LF(m,xind)=JIG_SamSize(m)*JIG_p_hat(m,xind);
         PS_LF(m,xind)=PS_SamSize(m)*PS_p_hat(m,xind);
         EFF_JIG_lengthfrq(m-indexMinmoLD+1,xind)=round(EFF_SamSize*(JIG_lengthfrq(m-indexMinmoLD+1,xind)/JIG_SamSize(m)));
         EFF_PS_lengthfrq(m-indexMinmoLD+1,xind)=round(EFF_SamSize*(PS_lengthfrq(m-indexMinmoLD+1,xind)/(PS_SamSize(m)+PS_delta(m))));
     };
  };
  cout<<"debug_10: "<<endl;

  JIG_logmult=0.0;
  PS_logmult=0.0;
  for(int i=indexMinmoLD;i<=nmos;i++) {
      if(EFF_Ssize_switch==0) {
            JIG_logmult+=factln(JIG_SamSize(i)); //log, factorial;
            PS_logmult+=factln(PS_SamSize(i))*PS_eta(i); //log, factorial;
      }
      else if(EFF_Ssize_switch==1) {
            JIG_logmult+=factln(EFF_SamSize); //log, factorial;
            PS_logmult+=factln(EFF_SamSize)*PS_eta(i); //log, factorial;
      };
      for(int xind=1;xind<=nlengths;xind++) { 
            if(EFF_Ssize_switch==0) {
                 JIG_logmult+=-1.0*factln((JIG_lengthfrq(i-indexMinmoLD+1,xind)))+(JIG_lengthfrq(i-indexMinmoLD+1,xind))*log(JIG_p_hat(i,xind));
                 PS_logmult+=-1.0*factln((PS_lengthfrq(i-indexMinmoLD+1,xind)))*PS_eta(i)+(PS_lengthfrq(i-indexMinmoLD+1,xind))*log(PS_p_hat(i,xind));
            }
            else if(EFF_Ssize_switch==1) {
                 JIG_logmult+=-1.0*factln((EFF_JIG_lengthfrq(i-indexMinmoLD+1,xind)))+(EFF_JIG_lengthfrq(i-indexMinmoLD+1,xind))*log(JIG_Catch(i,xind)/(sum(JIG_Catch(i))));
                 PS_logmult+=-1.0*factln((EFF_PS_lengthfrq(i-indexMinmoLD+1,xind)))*PS_eta(i)+(EFF_PS_lengthfrq(i-indexMinmoLD+1,xind))*log((PS_Catch(i,xind)+PS_delta(i))/((sum(PS_Catch(i)+PS_delta(i)))));
            };
      };
  };
  
  Q+=lambda1*(-1.0*(JIG_logmult)); //weight*(the negative multinomial likelihood);
  Q+=lambda2*(-1.0*(PS_logmult)); //weight*(the negative multinomial likelihood);
  cout<<"Q_1: "<<Q<<endl;
  maxloglike=0.0; 
  maxloglike+=JIG_logmult+PS_logmult;
  
  //part 2 of the objective function:log normal
  JIG_elem_obj2 = square(log(elem_div(JIG_yield,JIG_Yieldhat/1000.0))+JIG_sig2_yield/(2.0*square(JIG_Yieldhat/1000.0)))(indexMinmoLD,nmos);  //in MT
  PS_elem_obj2=square(log(elem_div(PS_yield+PS_delta,(PS_Yieldhat/1000.0)+PS_delta))+elem_prod(PS_eta,PS_sig2_yield/(2.0*square((PS_Yieldhat/1000.0)+PS_delta))))(indexMinmoLD,nmos);  //in MT
  
  
  Q+=(0.5*JIG_elem_obj2.size()*log(JIG_sig2_yield)+sum(JIG_elem_obj2)/(2.0*JIG_sig2_yield));
  Q+=(0.5*(PS_elem_obj2.size()-sum(PS_delta))*log(PS_sig2_yield)+sum(PS_elem_obj2)/(2.0*PS_sig2_yield));
  lognormal=0.0; 
  lognormal+=(-0.5*JIG_elem_obj2.size()*log(2*M_PI)-0.5*JIG_elem_obj2.size()*log(JIG_sig2_yield)-sum(JIG_elem_obj2)/(2.0*JIG_sig2_yield)); //M_PI=3.14192;
  lognormal+=(-0.5*(PS_elem_obj2.size()-sum(PS_delta))*log(2*M_PI)-0.5*(PS_elem_obj2.size()-sum(PS_delta))*log(PS_sig2_yield)-sum(PS_elem_obj2)/(2.0*PS_sig2_yield));
  cout<<"Q_2: "<<Q<<endl;
  maxloglike+=lognormal;
  
  if(LW_switch == 1) {
      //Part 3 for the allometric Length-weight relationship;
      // Weight|length ~ normal
      E_W=aWL*pow(lengthAL,bWL); 
      sig2AL=norm2(weightAL-E_W)/nAL;  
      Q+=lambda3*0.5*nAL*(1+log(sig2AL)); //negative loglikelihood
      
      maxloglike+=(-0.5*nAL*log(2*M_PI)-0.5*nAL*(1+log(sig2AL)));
      cout<<"Q_3: "<<Q<<endl;
  };
  
  if(SR_switch == 1 && choiceSR ==0 && Fecundity_switch == 1) {
      //Part 4 for the Fecundity;
      // Fecundity|length ~ lognormal
      E_FC=aFC*pow(lengthFC,bFC); 
      sig2FC=norm2(log(EggFC)-log(E_FC))/nFC;  
      Q+=lambda4*0.5*nFC*(1.0+log(sig2FC)); //negative loglikelihood
      
      maxloglike+=-0.5*nFC*log(2*M_PI)-0.5*nFC*log(sig2FC)-norm2(log(EggFC)-log(E_FC))/(2.0*sig2FC);
      cout<<"Q_4: "<<Q<<endl;
  };  
  
  /// AIC
  aic=-2.0*maxloglike+2.0*(nPARAMETERS); 
  
REPORT_SECTION
  report<<"Check the Results"<<endl;
  if(SR_switch == 0)
     report<<"    This is 'Original version'    "<<endl;
  else if(SR_switch == 1) {
        if(choiceSR  == 0)
            report<<"    This is 'S-R model using Fecundity version'    "<<endl;
        else if(choiceSR  == 1)
            report<<"    This is 'Beverton-Holt S-R model version'    "<<endl;
        else if(choiceSR  == 2)
            report<<"    This is 'Ricker S-R model version'    "<<endl;
  };
  report<<"                  "<<endl;
  report<<"Yr Month Recruits JIG_F_mo PS_F_mo JIG_Yield(Ton) JIG_Yieldhat(Ton) PS_Yield(Ton) PS_Yieldhat(Ton)"<<endl; 
  for(int i=1;i<=nmos;i++)  
     report<<YEAR(i)<<" "<<MONTH(i)<<" "<<Recruits(i)<<" "<<JIG_F_mo(i)<<" "<<PS_F_mo(i)<<" "<<JIG_yield(i)<<" "<<JIG_Yieldhat(i)/1000<<" "<<PS_yield(i)<<" "<<PS_Yieldhat(i)/1000<<endl; // Yield unit = ton
  if(SR_switch == 1 && Plus_switch == 1)
     report<<"2019"<<" "<<"1"<<" "<<Recruits(nmos+nProjection)<<" "<<"NA"<<" "<<"NA"<<" "<<"NA"<<" "<<"NA"<<" "<<"NA"<<" "<<"NA"<<endl; // Yield unit = ton
  report<<"         "<<endl;
  report<<"YEAR MONTH JIG_Effort PS_Effort Pop B.in.kg"<<endl;   
  for(int i=1;i<=nmos;i++)     
     report<<YEAR(i)<<" "<<MONTH(i)<<" "<<JIG_effort(i)<<" "<<PS_effort(i)<<" "<<Pop(i)<<" "<<B(i)<<" "<<endl;
  if(SR_switch == 1 && Plus_switch == 1)
          report<<"2019"<<" "<<"1"<<" "<<"NA"<<" "<<"NA"<<" "<<Pop(nmos+nProjection)<<" "<<B(nmos+nProjection)<<" "<<endl;
  report<<"         "<<endl;
  report<<"JIG_lengthfrq: "<<endl;
  report<<JIG_lengthfrq<<endl; 
  report<<"JIG_LF(): "<<endl;
  report<<JIG_LF<<endl;
  report<<"PS_lengthfrq: "<<endl;
  report<<PS_lengthfrq<<endl; 
  report<<"PS_LF(): "<<endl;
  report<<PS_LF<<endl;
  report<<"         "<<endl;
  if(LW_switch ==1) {
      report<<"weight given length (gram): "<<endl;
      report<<"Length Data: "<<lengthAL<<endl;
      report<<"Weight Data: "<<weightAL<<endl;
      report<<"predicted value: "<<E_W<<endl;
  };
  if(SR_switch ==1 && choiceSR == 0 && Fecundity_switch ==1) {
      report<<"Fecundity given length (No.): "<<endl;
      report<<"Length Data: "<<lengthFC<<endl;
      report<<"Fecundity Data: "<<EggFC<<endl;
      report<<"predicted value: "<<E_FC<<endl;
  };
  report<<"         "<<endl;
  report<<"The number of free parameters : "<<endl; 
  report<<nPARAMETERS<<endl;
  report<<"maximum gradiant : "<<endl; 
  report<<objective_function_value::gmax<<endl;
  report<<"Objective_function_value: "<<endl;
  report<<Q<<endl;
  report<<"AIC: "<<endl;
  report<<aic<<endl; 
  report<<"         "<<endl;
  report<<"Sample size for exploited  length composition data ("<<YEAR(indexMinmoLD)<<" yr "<<MONTH(indexMinmoLD)<<" bi-mo "<<" ~ "<<YEAR(nmos)<<" yr "<<MONTH(nmos)<<" bi-mo)"<<endl;
  report<<"JIG: "<<JIG_SamSize<<endl;
  report<<"PS: "<<PS_SamSize<<endl;
  report<<"         "<<endl;
  if(SR_switch == 1 && choiceSR == 0) {
      report<<"Time step || "<<"Eggs (No.)"<<endl;
      for(int m=indexMinmoLD; m<=nmos+nProjection;m++) {
            if(m < nmos+nProjection)
                   report<<YEAR(m)<<"yr "<<MONTH(m)<<"bi-mo || "<<Eggs(m)<<endl;
            else if(m == nmos+nProjection)
                   report<<"2019"<<"yr "<<"1"<<"bi-mo || "<<Eggs(m)<<endl;
      };
  };
  if(Plus_switch == 1)
    report<<"Time step || "<<"N: Age1 Age2 Age3 Age4 Age5 Age6+"<<endl;
  else if(Plus_switch == 0)
    report<<"Time step || "<<"N: Age1 Age2 Age3 Age4 Age5 Age6"<<endl;
  for(int m=1; m<=nmos+nProjection;m++) {
    if(SR_switch == 1) {
        if(m < nmos+nProjection)
            report<<YEAR(m)<<"yr "<<MONTH(m)<<"bi-mo || "<<trans(N)(m)<<endl;
        else if(m == nmos+nProjection)
            report<<"2019"<<"yr "<<"1"<<"bi-mo || "<<trans(N)(m)<<endl;
    }
    else if(SR_switch == 0) {
        report<<YEAR(m)<<"yr "<<MONTH(m)<<"bi-mo || "<<trans(N)(m)<<endl;
    };
  };
  report<<"         "<<endl;
  report<<"Time step || "<<"N_plus (No.)"<<endl;
  if(Plus_switch == 1) 
      for(int m=indexMinmoLD+1; m<=nmos+nProjection;m++) {
          if(SR_switch == 1) {
              if(m < nmos+nProjection)
                 report<<YEAR(m)<<"yr "<<MONTH(m)<<"bi-mo || "<<N_plus(m)<<endl;
              else if(m == nmos+nProjection)
                 report<<"2019"<<"yr "<<"1"<<"bi-mo || "<<N_plus(m)<<endl;
          }
          else if(SR_switch == 0)
              report<<YEAR(m)<<"yr "<<MONTH(m)<<"bi-mo || "<<N_plus(m)<<endl;
      };
  report<<"         "<<endl;
  if(Plus_switch == 1 && SR_switch == 1)
    report<<"Time step || "<<"Spawners: Age1 Age2 Age3 Age4 Age5 Age6+"<<endl;
  else if(Plus_switch == 0 && SR_switch == 1)
    report<<"Time step || "<<"Spawners: Age1 Age2 Age3 Age4 Age5 Age6"<<endl;
  if(SR_switch == 1) 
      for(int m=indexMinmoLD+1; m<=nmos+nProjection;m++) {
        if(m < nmos+nProjection)
            report<<YEAR(m)<<"yr "<<MONTH(m)<<"bi-mo || "<<trans(Spawners)(m)<<endl;
        if(m == nmos+nProjection)
            report<<"2019"<<"yr "<<"1"<<"bi-mo || "<<trans(Spawners)(m)<<endl;
      };
  report<<"         "<<endl;
  report<<"mu_r (cm): "<<mu_r<<endl;
  if(Plus_switch == 0) 
      report<<"Age: Age1 Age2 Age3 Age4 Age5 Age6"<<endl;
  if(Plus_switch == 1) 
      report<<"Age: Age1 Age2 Age3 Age4 Age5 Age6+"<<endl;
  report<<"SS[Age]: "<<SS<<endl;
  report<<"sigmaL (cm): "<<sigmaL<<endl; 
  report<<"         "<<endl;
  report<<"Length (cm): "<<x<<endl;   
  report<<"Mu (Length): "<<Mu<<endl;
  if(LW_switch == 1) 
      report<<"Wt (Length)_in_ kg: "<<Wt<<endl;
  if(Fecundity_switch == 1 && SR_switch ==1 && choiceSR == 0) 
      report<<"Fecundity (Length)_in_number: "<<fecundity<<endl;
  report<<"         "<<endl;
  if(SR_switch == 1 && Maturation_switch == 1) {
      report<<"Maturation(Length): "<<maturation<<endl;
      report<<"Spawning_M_at_Age 5 (Length): "<<-log(1.0-maturation*Sp_rate)<<endl;
      if(Plus_switch == 0)
          report<<"Spawning_M_at_Age 6 (Length): "<<-log(1.0-maturation)<<endl;
      if(Plus_switch == 1)
          report<<"Spawning_M_at_Age 6+ (Length): "<<-log(1.0-maturation)<<endl;
  };
  report<<"M(Length): "<<M<<endl; 
  if(SR_switch == 1 && Maturation_switch == 1) {
      report<<"Total_M_at_Age 1~4 (Length): "<<M<<endl; 
      report<<"Total_M_at_Age 5 (Length): "<<M-log(1.0-maturation*Sp_rate)<<endl;
      if(Plus_switch == 0)
          report<<"Total_M_at_Age 6 (Length): "<<M-log(1.0-maturation*Sp_rate)-log(1.0-maturation)<<endl; 
      if(Plus_switch == 1)
          report<<"Total_M_at_Age 6+ (Length): "<<M-log(1.0-maturation)<<endl; 
  };
  report<<"         "<<endl;
  report<<"JIG_Sel(Length): "<<JIG_Sel<<endl;
  report<<"PS_Sel(Length): "<<PS_Sel<<endl;
  report<<"         "<<endl;
  
  for(int c=1;c<=ncohorts+nProjection;c++) {
      if(c == 1)
      report<<"====== Start ======"<<endl;
      report<<"cohort: "<<c<<endl;
      report<<"f(Age, Length): "<<endl;
      report<<"               Length(cm): "<<x<<endl;
      for(int m=c;m<=c+indexMinmoLD-1;m++) {
         if(SR_switch == 1) {
            if(m < ncohorts+nProjection) {
                if(Plus_switch == 0)
                    report<<"   "<<YEAR(m)<<"yr "<<MONTH(m)<<"bi-mo ||"<<" Age"<<(m-c+1)<<": "<<NL(m-(c-1),m)/N(m-(c-1),m)<<endl;
                else if(Plus_switch == 1 && (m-c+1) != nages) 
                    report<<"   "<<YEAR(m)<<"yr "<<MONTH(m)<<"bi-mo ||"<<" Age"<<(m-c+1)<<": "<<NL(m-(c-1),m)/N(m-(c-1),m)<<endl;
                else if(Plus_switch == 1 && (m-c+1) == nages) 
                    report<<"   "<<YEAR(m)<<"yr "<<MONTH(m)<<"bi-mo ||"<<" Age"<<(m-c+1)<<"+: "<<NL(m-(c-1),m)/N(m-(c-1),m)<<endl;
            }
            else if(m == ncohorts+nProjection) {
                if(Plus_switch == 0)
                    report<<"   "<<"2019"<<"yr "<<"1"<<"bi-mo ||"<<" Age"<<(m-c+1)<<": "<<NL(m-(c-1),m)/N(m-(c-1),m)<<endl;
                else if(Plus_switch == 1 && (m-c+1) != nages) 
                    report<<"   "<<"2019"<<"yr "<<"1"<<"bi-mo ||"<<" Age"<<(m-c+1)<<": "<<NL(m-(c-1),m)/N(m-(c-1),m)<<endl;
                else if(Plus_switch == 1 && (m-c+1) == nages) 
                    report<<"   "<<"2019"<<"yr "<<"1"<<"bi-mo ||"<<" Age"<<(m-c+1)<<"+: "<<NL(m-(c-1),m)/N(m-(c-1),m)<<endl;
             };
         }
         else if(SR_switch == 0) 
            if(m <= ncohorts+nProjection) {
                if(Plus_switch == 0)
                    report<<"   "<<YEAR(m)<<"yr "<<MONTH(m)<<"bi-mo ||"<<" Age"<<(m-c+1)<<": "<<NL(m-(c-1),m)/N(m-(c-1),m)<<endl;
                else if(Plus_switch == 1 && (m-c+1) != nages) 
                    report<<"   "<<YEAR(m)<<"yr "<<MONTH(m)<<"bi-mo ||"<<" Age"<<(m-c+1)<<": "<<NL(m-(c-1),m)/N(m-(c-1),m)<<endl;
                else if(Plus_switch == 1 && (m-c+1) == nages) 
                    report<<"   "<<YEAR(m)<<"yr "<<MONTH(m)<<"bi-mo ||"<<" Age"<<(m-c+1)<<"+: "<<NL(m-(c-1),m)/N(m-(c-1),m)<<endl;
            }; 
       };
       for(int m=c;m<=c+indexMinmoLD-1;m++)
           if(m <= ncohorts+nProjection) {
                if(Plus_switch == 0)
                    report<<"Expected length of Age ["<<(m-c+1)<<"]: "<<sum(elem_prod(x,NL(m-(c-1),m)/N(m-(c-1),m)))<<" cm"<<endl;
                else if(Plus_switch == 1 && (m-c+1) != nages) 
                    report<<"Expected length of Age ["<<(m-c+1)<<"]: "<<sum(elem_prod(x,NL(m-(c-1),m)/N(m-(c-1),m)))<<" cm"<<endl;
                else if(Plus_switch == 1 && (m-c+1) == nages) 
                    report<<"Expected length of Age ["<<(m-c+1)<<"+]: "<<sum(elem_prod(x,NL(m-(c-1),m)/N(m-(c-1),m)))<<" cm"<<endl;
           };
      report<<"NL(Age, Cohort, Length): "<<endl;
      report<<"           Length(cm): "<<x<<endl;
      for(int m=c;m<=c+indexMinmoLD-1;m++) {
         if(SR_switch == 1) {
            if(m < ncohorts+nProjection) {
                if(Plus_switch == 0)
                    report<<"   "<<YEAR(m)<<"yr "<<MONTH(m)<<"bi-mo ||"<<" Age"<<(m-c+1)<<": "<<NL(m-(c-1),m)<<endl;
                else if(Plus_switch == 1 && (m-c+1) != nages) 
                    report<<"   "<<YEAR(m)<<"yr "<<MONTH(m)<<"bi-mo ||"<<" Age"<<(m-c+1)<<": "<<NL(m-(c-1),m)<<endl;
                else if(Plus_switch == 1 && (m-c+1) == nages) 
                    report<<"   "<<YEAR(m)<<"yr "<<MONTH(m)<<"bi-mo ||"<<" Age"<<(m-c+1)<<"+: "<<NL(m-(c-1),m)<<endl;
            }
            else if(m == ncohorts+nProjection) {
                if(Plus_switch == 0)
                    report<<"   "<<"2019"<<"yr "<<"1"<<"bi-mo ||"<<" Age"<<(m-c+1)<<": "<<NL(m-(c-1),m)<<endl;
                else if(Plus_switch == 1 && (m-c+1) != nages) 
                    report<<"   "<<"2019"<<"yr "<<"1"<<"bi-mo ||"<<" Age"<<(m-c+1)<<": "<<NL(m-(c-1),m)<<endl;
                else if(Plus_switch == 1 && (m-c+1) == nages) 
                    report<<"   "<<"2019"<<"yr "<<"1"<<"bi-mo ||"<<" Age"<<(m-c+1)<<"+: "<<NL(m-(c-1),m)<<endl;
             };
         }
         else if(SR_switch == 0) 
            if(m <= ncohorts+nProjection) {
                if(Plus_switch == 0)
                    report<<"   "<<YEAR(m)<<"yr "<<MONTH(m)<<"bi-mo ||"<<" Age"<<(m-c+1)<<": "<<NL(m-(c-1),m)<<endl;
                else if(Plus_switch == 1 && (m-c+1) != nages) 
                    report<<"   "<<YEAR(m)<<"yr "<<MONTH(m)<<"bi-mo ||"<<" Age"<<(m-c+1)<<": "<<NL(m-(c-1),m)<<endl;
                else if(Plus_switch == 1 && (m-c+1) == nages) 
                    report<<"   "<<YEAR(m)<<"yr "<<MONTH(m)<<"bi-mo ||"<<" Age"<<(m-c+1)<<"+: "<<NL(m-(c-1),m)<<endl;
            }; 
       };
       if(c < ncohorts+nProjection)
           report<<"============>>> Next cohort"<<endl;
       else if(c == ncohorts+nProjection)
           report<<"====== Finish ======"<<endl;
  };
  
RUNTIME_SECTION
  maximum_function_evaluations 100,150,300,1000
  convergence_criteria .001,.0001,1e-7
    
TOP_OF_MAIN_SECTION
  gradient_structure::set_MAX_NVAR_OFFSET(1000);  //maximum number of depdendent variables of 400 exceeded 
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(1000);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(100000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(1000000);
  arrmblsize=900000;
  
GLOBALS_SECTION
  #include <admodel.h>
  #include <math.h>
  #include <stdio.h>
  #include <stddef.h>
  #include <stdlib.h>

