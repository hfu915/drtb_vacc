%% Setup_model_ct1.m
% set up and model input parameters
% update on 01/12/2020

% Select country for evaluation
% icty = 1;

% load country data
DataDmg = readtable('dataHBC_dmg.xlsx'); % units in thousands
DataTbb = readtable('dataHBC_tbb.xlsx');
DataTbbYrs = readtable('dataHBC_tbb_yrs.xlsx');
DataCtt = readtable('dataHBC_ctt.xlsx');

%% set up compartments
gps.ages    = {'ch','ad'};
gps.vaccs   = {'v0','v1','v2'};                                           % unvaccinated, vaccinated and immuned, vaccinated but waned
gps.strains = {'DS','DR'};
gps.sectors = {'pu','pr'};

states0 = {'U'};                                                           % Structured only by age, vacc
states1 = {'L0','L1','I','E','Rlo','Rhi','R'};                             % + by strain
states2 = {'Dx','Tx','Tx2'};                                               % + by sector

[i, s, d, lim] = get_addresses({states0, gps.ages, gps.vaccs}, [], [], [], 0);
[i, s, d, lim] = get_addresses({states1, gps.ages, gps.vaccs, gps.strains}, i, s, d, lim);
[i, s, d, lim] = get_addresses({states2, gps.ages, gps.vaccs, gps.strains, gps.sectors}, i, s, d, lim);
d = char(d);

% Include the auxiliaries
names = {'inc','mort'};
lgths = [    6,     3];                                                    % Total, ch DS, ad DS, ch MDR, ad MDR
for ii = 1:length(names)
    inds = lim+[1:lgths(ii)];
    i.aux.(names{ii}) = inds;
    lim = inds(end);
end

i.aux.ntfsec = i.nstates + sum(lgths) + [1 2]; 
i.aux.newinf = i.nstates + sum(lgths) + [3 4];
i.allcomp = i.nstates + sum(lgths) + 4;

% Specify matrices for selection
% New infection: 1.DS, 2.DR    
tmp = zeros(2,i.nstates);
tmp(1,intersect(s.L0,s.DS)) = 1;  
tmp(2,intersect(s.L0,s.DR)) = 1;  
i.agg.newinf = sparse(tmp);
    
tmp = zeros(i.nstates);
tmp(s.L0,:) = 1;
i.sel.newinf = sparse(tmp - diag(diag(tmp)));

s.infectious = [s.I, s.Dx, s.E, intersect(s.Tx, s.DR)];
s.prev_TB    = unique([s.I, s.Dx, s.E, s.Tx, s.Tx2]);


%% Specify parameters 
% Calibration and evaluation of uncertainty
r.birth        = -log(1-DataDmg.cbirth(icty));                  %x(1);     % Birth rate 
r.ageing       = 1/15;                                          %x(2);     % Ageing rate
r.beta         = [15 8];                                        %x(3:4);   % Infection rate for DS- and DR-TB                              
betaRt_prm     = 0.98;                                          %x(5);     % Slope for declining r.beta
r.careseeking  = 1;                                             %x(6);     % Initial health-care seeking rate
Tx_init_pu     = DataTbb.txcov(icty)/100;                       %x(7);     % Average proportion of treatment initiation, WHO CDR
reactiv        = 0.001;                                         %x(8);     % TB reactivation rate in adults
breakdwn       = 0.0527;                                        %x(9);     % Primary TB progression rate in adults  
pgrCh          = 0.5;                                           %x(10);    % Reduced risk of TB progression in children comapred to adults
p.red_sus      = 0.5;                                           %x(11);    % Reduced risk of TB infection in LTBI
p.red_infectCh = 0.824;                                         %x(12);    % Reduced susceptability in children comapred to adults
r.mort_TB      = 1/6;                                           %x(13);    % TB mortality rate
p.pu           = 0.5;                                           %x(14);    % Current proportion of care-seeking in public sectors  
r.DR_acqu      = 0.01;                                          %x(15);    % Rate of acquiring resistance during treatment
p.inferior_pr  = 0.6;                                           %x(16);    % Reduced quality of care in private sectors compared to public ones
r.careseeking2 = 12;                                            %x(17);    % Secondary health-care seeking rate
r.relapse      = [0.032 0.14 0.0015];                           %x(18:20); % Recurrence rates for success treatment, defaulters, stabalised state after to years

% Scale vector to improving mixing for calibration
p.AdjScale = [100 100 1/10 1 10 1 10 1e3 100 .....
    10 10 10 10 10 100 10 1 100 10 1e3];

% Natural history
p.betaRt_yr     = 1970;                                                    % Beginning year of declining r.beta
r.reactivation  = reactiv*[1; 1]*[1 1 1];                                  % Reactivation rate by age and vaccine status
r.breakdown     = breakdwn*[1; 1]*[1 1 1];                                 % Primary progression rate by age and vaccine status
r.stabilisation = 1/2; 

mcontact  = [DataCtt.ch_ch(icty) DataCtt.ch_ad(icty);....                  
             DataCtt.ad_ch(icty) DataCtt.ad_ad(icty)];                     % Contact rates (susceptible, infected)
p.contact = mcontact/sum(sum(mcontact))*2;                                 % Standardised contact rates (make beta comparable with the homogeneous setting)

r.mort          = [-log(1-DataDmg.Ch_death(icty)/DataDmg.Ch_pop(icty)) ... % Natural-cause mortality (children)
                   1/DataDmg.lifeexpAd(icty)];                             % Natural-cause mortality (adults)
r.self_cure     = r.mort_TB;                                               % Rate of self recovery (assumed to be r.mort_TB)

% Diagnosis and linkage to treatment (sector-specific)
r.Dx           = 52;
p.Tx_initFL    = Tx_init_pu*[1 p.inferior_pr];
p.Tx_initSL    = Tx_init_pu*[1 0];
p.DR_rec     = DataTbb.pct_rrtest_wt(icty)*[1 0];                          % Proportion of knowing DR state at initial health seeking
p.SL_trans   = [0.85 0];                                                   % Proportion of DR cases switching to SL after FL treatment, by sector

% FL treatment outcomes
p.FLlost     = DataTbb.pct_FL_lost(icty);                                  % Proportion of loss-to-follow-up during FL treatment
p.FLdied     = DataTbb.pct_FL_died(icty);                                  % Proportion of death during FL treatment
p.FLfail     = DataTbb.pct_FL_fail(icty);                                  % Proportion of failure during FL treatment
p_FLlostPuPr = p.FLlost*[1 2-p.inferior_pr];                               % Proportion of loss-to-follow-up during FL treatment, private sector
p_FLsuccPuPr = 1 - p.FLdied - p.FLfail - p_FLlostPuPr;                     % Proportion of success during FL treatment
r.txFLcom    = 2;                                                          % 1/Duration of FL treatment completion
r_txFLnotcom = r.txFLcom./(p_FLsuccPuPr+p.FLfail)-r.txFLcom;               % Rate of FL treatment completion
r.txFLdef    = p_FLlostPuPr.*(r.txFLcom+r_txFLnotcom);                     % Rate of loss-to-follow-up during FL treatment 
r.txFLdie    = p.FLdied.*(r.txFLcom+r_txFLnotcom);                         % Rate of death during FL treatment
p.cureFL     = (p_FLsuccPuPr./(p_FLsuccPuPr+p.FLfail));                    % Proportion of success among those complete FL treatment

% SL treatment outcomes
p_SLlost     = DataTbb.pct_SL_lost(icty);                                  % Proportion of loss-to-follow-up during SL treatment
p_SLdied     = DataTbb.pct_SL_died(icty);                                  % Proportion of death during SL treatment
p_SLsucc     = DataTbb.pct_SL_succ(icty);                                  % Proportion of success during SL treatment
p_SLfail     = DataTbb.pct_SL_fail(icty);                                  % Proportion of failure during SL treatment
r.txSLcom    = 0.5;                                                        % 1/Duration of SL treatment completion
r_txSLnotcom = r.txSLcom/(p_SLsucc+p_SLfail)-r.txSLcom ;                   % Rate of SL treatment completion
r.txSLdef    = p_SLlost*(r.txSLcom+r_txSLnotcom);                          % Rate of loss-to-follow-up during SL treatment 
r.txSLdie    = p_SLdied*(r.txSLcom+r_txSLnotcom);                          % Rate of death during SL treatment
p.cureSL     = (p_SLsucc/(p_SLsucc+p_SLfail))*[1 0];                       % Proportion of success among those complete SL treatment, by sector

% Interventions
p.VE            = 0;                                                       % Vaccine efficacy in reducing susceptibility
r.vacc_catchup  = [0 0];                                                   % Coverage amongst newborns
p.vacc_birthday = [0 0];                                                   % Coverage amongst 15yo: LTBI vs whole population
r.waning        = 0.1;


prm.p = p; prm.r = r;
ref.i = i; ref.s = s; ref.d = d;

%% Load data with 95%CIs
% Demographics
CtyDmg = DataDmg(icty,:);
data.pop1970 = CtyDmg.yr1970; % 130148.65 in thousands
pop2018 = CtyDmg.yr2018; % 145530.09 in thousands
data.popgrowth = (pop2018/data.pop1970)*[0.8 1 1.2];  % 1.1182
data.pctCh = (CtyDmg.Ch_pop/pop2018)*[0.8 1 1.2]; % 0.1763

% TB burden
% Single timepoints
CtyTbb = DataTbb(icty,:);
data.incdALch = [CtyTbb.incrtCh_lo    CtyTbb.incrtCh     CtyTbb.incrtCh_hi ]; 
data.incdALad = [CtyTbb.incrtAd_lo    CtyTbb.incrtAd     CtyTbb.incrtAd_hi ];  
data.incdDR   = [CtyTbb.incrtDR_lo    CtyTbb.incrtDR     CtyTbb.incrtDR_hi ];  
data.ntfPu    = CtyTbb.ntfrtPub*[0.8 1 1.2];                                       
% Multiple timepoints
CtyTbbYrs     = DataTbbYrs(string(DataTbbYrs.iso3)==CtyISO3{icty},:);
data.incdAL   = [CtyTbbYrs.e_inc_100k_lo   CtyTbbYrs.e_inc_100k    CtyTbbYrs.e_inc_100k_hi]; 
data.mortAL   = [CtyTbbYrs.e_mort_100k_lo   CtyTbbYrs.e_mort_100k   CtyTbbYrs.e_mort_100k_hi]; 
