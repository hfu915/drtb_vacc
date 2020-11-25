%% Get_vacc_results_ct1.m
% get projections of vaccination and intervention scenarios
% update on 12/11/2020

% --- vaccination strategies 
% (1) Regular vaccination: all 15-year-old at their birthdays, 2020-2035
% (2) Campaign vaccination: vaccinate all adults, 2025-2026, 2030-2031

% --- vaccination scenarios
% interventions: A - improved RR-TB care, B - M72-like vaccination
% (1) A - no,  B - no  (status quo)
% (2) A - yes, B - no
% (3) A - no,  B - yes
% (4) A - yes, B - yes

clear

CtyISO3 = {'IND' 'CHN' 'RUS' 'PAK' 'IDN' ....
 'NGA' 'PHL' 'UKR' 'MMR' 'ZAF' ....
 'VNM' 'MOZ' 'COD' 'BGD' 'PRK' ....
 'KAZ' 'UZB' 'SOM' 'THA' 'AGO' ....
 'PER' 'KGZ' 'KEN' 'PNG' 'TJK' ....
 'ETH' 'ZWE' 'BLR' 'MDA' 'AZE'};

% select category 1 countries (#7)
CtyCt1 = [1 4 5 7 9 14 26];
nCt1 = length(CtyCt1);

% load data for vaccination coverages
DataVacCov = readtable('dataHBC_VacCov.xlsx');

% select a level of vaccine efficacy for evaluation
sel_ve = 0.5;
evaVE = num2str(sel_ve*100); 

for iord = 1:1 % running projectionos for India only %1:nCt1 
    icty = CtyCt1(iord);
    
    % load dataset
    evaCty = CtyISO3{icty};
    load(['MCMC_posterior_' evaCty '.mat']);
    nPostS = size(xs,1);

    % assign coverages for two vaccination programmes
    secatt_cov = DataVacCov.secattcov(icty);
    camp_cov = DataVacCov.campcov(icty);
    
    % add a counter for vaccination doses
    i.allcomp = i.allcomp + 1;
    i.aux.vdose = i.allcomp;

    tmp = zeros(i.nstates);
    tmp(s.v1,intersect(s.v0,[s.L0 s.L1])) = 1;
    i.sel.vdose = sparse(tmp - diag(diag(tmp)));
    
    % Intervention A: set up parameters for DST improvement (access=85%)
    a_pDR_rec  = max(0.85, p.DR_rec(1));

    % Intervention A: set up parameters for and SL improvement (success=75%) in public sectors
    a_p_SLsucc   = max(0.75, p_SLsucc);                                                       
    gap_SLsucc   = a_p_SLsucc - p_SLsucc;
    a_p_SLlost   = p_SLlost - gap_SLsucc*p_SLlost/(1-p_SLsucc);
    a_p_SLdied   = p_SLdied - gap_SLsucc*p_SLdied/(1-p_SLsucc);
    a_p_SLfail   = p_SLfail - gap_SLsucc*p_SLfail/(1-p_SLsucc);
    a_r_txSLnotcom = r.txSLcom/(a_p_SLsucc+a_p_SLfail)-r.txSLcom ;                   

    Allsoln = zeros(length(2019:2035)+1, i.allcomp, 4, nPostS);                               % dim: (3) vacc                      
    mk = round(nPostS/20);

    for ix = 1:nPostS
        if mod(ix,mk)==0; fprintf('%2.1f%% done\n',100*(ix/nPostS)); end

        x = xs(ix,:)./p.AdjScale;

        r.birth        = x(1);  
        r.ageing       = x(2);
        r.beta         = x(3:4);
        betaRt_prm     = x(5);
        r.careseeking  = x(6);
        Tx_init_pu     = x(7);
        reactiv        = x(8);  
        breakdwn       = x(9);  
        pgrCh          = x(10); 
        p.red_sus      = x(11);  
        p.red_infectCh = x(12); 
        r.mort_TB      = x(13); 
        p.pu           = x(14); 
        r.DR_acqu      = x(15); 
        p.inferior_pr  = x(16);
        r.careseeking2 = x(17);  
        r.relapse      = x(18:20);
        
        r.self_cure  = r.mort_TB;

        r.betaRt_fun    = @(t)betaRt_prm^(t-p.betaRt_yr);
        r.reactivation  = reactiv*[pgrCh; 1]*[1 1 1];                          % dim: (1)age (2)vacc
        r.breakdown     = breakdwn*[pgrCh; 1]*[1 1 1];

        p.Tx_initFL     = Tx_init_pu*[1 p.inferior_pr];                        % dim: pu/pr
        p.Tx_initSL     = Tx_init_pu*[1 0];                                    % dim: pu/pr
        p_FLlostPuPr = p.FLlost*[1 2-p.inferior_pr];                           % Proportion of loss-to-follow-up during FL treatment, private sector
        p_FLsuccPuPr = 1 - p.FLdied - p.FLfail - p_FLlostPuPr;                 % Proportion of success during FL treatment                      
        r_txFLnotcom = r.txFLcom./(p_FLsuccPuPr+p.FLfail)-r.txFLcom;           % Rate of FL treatment completion
        r.txFLdef    = p_FLlostPuPr.*(r.txFLcom+r_txFLnotcom);                 % Rate of loss-to-follow-up during FL treatment 
        r.txFLdie    = p.FLdied.*(r.txFLcom+r_txFLnotcom);                     % Rate of death during FL treatment
        p.cureFL     = (p_FLsuccPuPr./(p_FLsuccPuPr+p.FLfail));                % Proportion of success among those complete FL treatment

        %% --- Set up parameters for intervention scenarios 
        p.equpop = 0;
        % vaccine scenario (1) - status quo
        p.vacc_birthday = [0 0];
        r.vacc_catchup  = [0 0];
        Mpjt = make_model_ct1(p, r, i, s, gps, []);

        % vaccine scenario (2) - improved RR-TB care only (A)
        ar = r; ap = p;
        ar.txSLdef    = a_p_SLlost*(r.txSLcom+a_r_txSLnotcom);                          
        ar.txSLdie    = a_p_SLdied*(r.txSLcom+a_r_txSLnotcom);                          
        ap.cureSL     = (a_p_SLsucc/(a_p_SLsucc+a_p_SLfail))*[1 0]; 
        ap.DR_rec(1)  = a_pDR_rec;
        Mtximp = make_model_ct1(ap, ar, i, s, gps, []);    

        % vaccine scenario (3) - M72-like vaccine only (B)
        br = r; bp = p;
        br.vacc_catchup  = [0 camp_cov];
        bp.vacc_birthday = [0 secatt_cov];

        br.reactivation  = reactiv*[pgrCh; 1]*[1 (1-sel_ve) 1];    % dim: (1)age (2)vacc
        br.breakdown     = breakdwn*[pgrCh; 1]*[1 (1-sel_ve) 1];
        bp.VE = 0;

        MvacUp = make_model_ct1(bp, br, i, s, gps, [s.L0, s.L1]);

        br.vacc_catchup  = [0 0];
        MvacDn = make_model_ct1(bp, br, i, s, gps, [s.L0, s.L1]);

        % vaccine scenario (4) - improved RR-TB care + M72-like vaccine (A+B)
        abr = ar; abp = ap;
        abr.vacc_catchup  = [0 camp_cov];
        abp.vacc_birthday = [0 secatt_cov];

        abr.reactivation  = reactiv*[pgrCh; 1]*[1 (1-sel_ve) 1];    % dim: (1)age (2)vacc
        abr.breakdown     = breakdwn*[pgrCh; 1]*[1 (1-sel_ve) 1];
        abp.VE = 0;

        MvacUp_tx = make_model_ct1(abp, abr, i, s, gps, [s.L0, s.L1]);

        abr.vacc_catchup  = [0 0];
        MvacDn_tx = make_model_ct1(abp, abr, i, s, gps, [s.L0, s.L1]);


        %% --- Incorporate parameters and initiate simulation

        % future projections - before introduction of vaccination
        init_pjt = [ResSolnEnd(ix,:) 0];          % initial distribution of population from previous results
        getsol_rn0 = @(Mtx,r0,p0) ode15s(@(t,in) goveqs_scaleup_int_ct1(t, in, Mpjt, Mtx, [2020 2023], i, s, r0, p0), 2019:1:2025, init_pjt, odeset('NonNegative',1:i.nstates));

        % first-round campaign
        getsol_rn1_up = @(Mup,Mdn,init_r1_up,r1,p1) ode15s(@(t,in) goveqs_scaleup_ct1(t, in, Mdn, Mup, [2025 2026], i, s, r1, p1), 2025:0.1:2026, init_r1_up, odeset('NonNegative',1:i.nstates));
        getsol_rn1_dn = @(Mup,Mdn,init_r1_dn,r1,p1) ode15s(@(t,in) goveqs_scaleup_ct1(t, in, Mup, Mdn, [2026 2027], i, s, r1, p1), 2026:0.1:2028, init_r1_dn, odeset('NonNegative',1:i.nstates));

        % second-round campaign 
        getsol_rn2_up = @(Mup,Mdn,init_r2_up,r2,p2) ode15s(@(t,in) goveqs_scaleup_ct1(t, in, Mdn, Mup, [2030 2031], i, s, r2, p2), 2028:0.1:2031, init_r2_up, odeset('NonNegative',1:i.nstates));
        getsol_rn2_dn = @(Mup,Mdn,init_r2_dn,r2,p2) ode15s(@(t,in) goveqs_scaleup_ct1(t, in, Mup, Mdn, [2031 2032], i, s, r2, p2), 2031:0.1:2036, init_r2_dn, odeset('NonNegative',1:i.nstates));

        % Simulate the models
        allMtximp = {Mpjt, Mtximp,   Mpjt,    Mtximp};
        allr_rn0  = {   r,     ar,      r,        ar};
        allp_rn0  = {   p,     ap,      p,        ap};
        allMvacUp = {Mpjt, Mtximp, MvacUp, MvacUp_tx};
        allMvacDn = {Mpjt, Mtximp, MvacDn, MvacDn_tx}; 
        allr_rn12 = {   r,     ar,     br,       abr};
        allp_rn12 = {   p,     ap,     bp,       abp};

        for im = 1:length(allMvacUp)
            [T0,auxTmp0] = getsol_rn0(allMtximp{im},allr_rn0{im},allp_rn0{im});
            [T1,auxTmp1] = getsol_rn1_up(allMvacUp{im},allMvacDn{im},auxTmp0(end,:),allr_rn12{im},allp_rn12{im});
            [T2,auxTmp2] = getsol_rn1_dn(allMvacUp{im},allMvacDn{im},auxTmp1(end,:),allr_rn12{im},allp_rn12{im});
            [T3,auxTmp3] = getsol_rn2_up(allMvacUp{im},allMvacDn{im},auxTmp2(end,:),allr_rn12{im},allp_rn12{im});
            [T4,auxTmp4] = getsol_rn2_dn(allMvacUp{im},allMvacDn{im},auxTmp3(end,:),allr_rn12{im},allp_rn12{im});

            % combine all results and intrapolate to get estimates of specific years
            Tall = [T0(1:end); T1(2:end); T2(2:end); T3(2:end); T4(2:end)];
            auxTmpAll = cat(1,auxTmp0(1:end,:),auxTmp1(2:end,:),auxTmp2(2:end,:),auxTmp3(2:end,:),auxTmp4(2:end,:)); 
            Allsoln(:,:,im,ix) = interp1(Tall,auxTmpAll,2019:2036); 

            clear auxTmp0 auxTmp1 auxTmp2 auxTmp3 auxTmp4 auxTmpAll Tall
        end
    end

    save(['Vacc_Tx_' CtyISO3{icty} '_ve' evaVE '.mat'],'Allsoln','s','i','data');
    fprintf('country finished: %s, VE: %s%%\n',evaCty,evaVE)
    clearvars -except CtyCt1 CtyISO3 sel_ve evaVE DataVacCov

end
