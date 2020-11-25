%% Get_summary_results.m
% get global estimates and summarise overall results of interventions for publication
% conduct after getting country-specific projections (Vacc_Tx_*_ve50.mat)

clear

% load data for antibiotic consumption analysis
DataCough = readtable('dataHBC_cough_IHME.xlsx');
DataAntib = readtable('dataHBC_antibiotic_DHS.xlsx');

CtyISO3 = {'IND' 'CHN' 'RUS' 'PAK' 'IDN' ....
 'NGA' 'PHL' 'UKR' 'MMR' 'ZAF' ....
 'VNM' 'MOZ' 'COD' 'BGD' 'PRK' ....
 'KAZ' 'UZB' 'SOM' 'THA' 'AGO' ....
 'PER' 'KGZ' 'KEN' 'PNG' 'TJK' ....
 'ETH' 'ZWE' 'BLR' 'MDA' 'AZE'};

Gestims    = zeros(30,4,200,6);    % Dims: 1.country 2.int 3.sample 4.indicators
Prcestims  = zeros(30,3,3,9);      % Dims: 1.country 2.int 3.percentile 4.indicators

% determine the scope of cough-related diseases
CoughR = {'ExURI','All','ExULRI','LRITB'};
ico = 1;

for icty = 1:1 % processing India only % 1:30
    
    load(['Vacc_Tx_' CtyISO3{icty} '_ve50.mat'],'Allsoln','i','s','data'); 
    Aallsoln = (1/2)*(Allsoln(1:(end-1),:,:,:) + Allsoln(2:end,:,:,:));    % mid-time estimates
    Dallsoln = diff(Allsoln,[],1);
    popn = data.pop1970*1e3; 

    selYrs = (2019:2035)-2019+1;
    Amat   = squeeze(sum(Aallsoln(selYrs,:,:,:),1));
    Dmat   = squeeze(sum(Dallsoln(selYrs,:,:,:),1));                       % Dims: 1.States 2.Intvn 3.Samples

    % (1) Averted TB-DR cases, proportion of aversion
    CumInc   = squeeze(sum(Dmat(i.aux.inc([4 5]),:,:),1));          % Dims: 1.Intvn 2.Samples
    NCumCase = CumInc*popn;
    Ncase    = (CumInc(1,:) - CumInc(2:end,:))*popn;
    NavtInc  = prctile(Ncase,[2.5,50,97.5],2);
    Pcase    = Ncase./(CumInc(1,:)*popn);
    PavtInc = prctile(Pcase,[2.5,50,97.5],2);

    % (2) Averted TB-DR deaths
    CumMor   = squeeze(Dmat(i.aux.mort(3),:,:));
    NCumDeath = CumMor*popn;
    Ndeath   = (CumMor(1,:) - CumMor(2:end,:))*popn;
    NavtMor = prctile(Ndeath,[2.5,50,97.5],2);
    Pdeath   = Ndeath./(CumMor(1,:)*popn);
    PavtMor = prctile(Pdeath,[2.5,50,97.5],2);

    % (3) Patient-months under FL treatment
    CumPMTx1  = squeeze(sum(Amat(s.Tx,:,:),1));
    NCumPMTx1 = 12*popn*CumPMTx1;
    NPMTx1    = 12*popn*(CumPMTx1(1,:) -  CumPMTx1(2:end,:));
    NavtTx1   = prctile(NPMTx1,[2.5,50,97.5],2);

    % (4) Patient-month under SL treatment
    CumPMTx2  = squeeze(sum(Amat(s.Tx2,:,:),1));
    NCumPMTx2 = 12*popn*CumPMTx2;
    NPMTx2    = 12*popn*(CumPMTx2(1,:) -  CumPMTx2(2:end,:));
    NavtTx2   = prctile(NPMTx2,[2.5,50,97.5],2);
    PPMTx2    = NPMTx2./(12*popn*CumPMTx2(1,:));
    PavtTx2   = prctile(PPMTx2,[2.5,50,97.5],2);

    % (5) Proportion of cough aetiology reduction
    pctCoughTB = DataCough.TBMed_Nm(icty)/DataCough.([CoughR{ico} 'Med_Nm'])(icty);
    pTB0 = pctCoughTB*(0.8 + 0.4*rand(1,size(Allsoln,4)));                 % +-20%, uniform
    pTB1 = pTB0.*Pcase;                                                    % Dims 1.Intvn 2.nPostP 
    RCE = prctile(pTB1, [2.5,50,97.5],2);

    % (6) Reduction in number of TB cases with antibiotic treatment
    pctAbiUse  = DataAntib.AbiUse(icty)/100;
    pAbi = pctAbiUse*(0.8 + 0.4*rand(1,size(Allsoln,4))); % +-20%
    RdNumAbiUse = Ncase.*pctAbiUse;
    RAB = prctile(RdNumAbiUse, [2.5,50,97.5], 2);

    % (7-10) Total TB cases, DR-TB cases, in 2020 and 2035
    selyr2020 = 2020-2019+1;  selyr2035 = 2035-2019+1;
    AlTBInc2020  = popn*squeeze(Dallsoln(selyr2020,i.aux.inc(1),:,:));     % Dims: 1.Intvn 2.Samples
    AlTBInc2035  = popn*squeeze(Dallsoln(selyr2035,i.aux.inc(1),:,:));    
    DrTBInc2020  = popn*squeeze(sum(Dallsoln(selyr2020,i.aux.inc([4 5]),:,:),2)); 
    DrTBInc2035  = popn*squeeze(sum(Dallsoln(selyr2035,i.aux.inc([4 5]),:,:),2)); 

    % Incremental effects of vaccine to RR-TB management on averted burden
    Ncase_Vx2Tx = CumInc(2,:) - CumInc(4,:);
    Pcase_Vx2Tx = Ncase_Vx2Tx./CumInc(1,:);
    PavtInc_Vx2Tx = prctile(Pcase_Vx2Tx,[2.5,50,97.5]);

    Ndeath_Vx2Tx = CumMor(2,:) - CumMor(4,:);
    Pdeath_Vx2Tx = Ndeath_Vx2Tx./CumMor(1,:);
    PavtMor_Vx2Tx = prctile(Pdeath_Vx2Tx,[2.5,50,97.5]);

    NPMTx2_Vx2Tx = CumPMTx2(2,:) - CumPMTx2(4,:);
    PPMTx2_Vx2Tx = NPMTx2_Vx2Tx./CumPMTx2(1,:);
    PavtTx_Vx2Tx = prctile(PPMTx2_Vx2Tx,[2.5,50,97.5]);

    Avt_Vx2Tx = [PavtInc_Vx2Tx; PavtMor_Vx2Tx; PavtTx_Vx2Tx];


    % get raw estimates of all posteiror samples
    Gestims(icty,:,:,1)   = NCumCase;
    Gestims(icty,:,:,2)   = NCumDeath;
    Gestims(icty,:,:,3)   = NCumPMTx1;
    Gestims(icty,:,:,4)   = NCumPMTx2;
    Gestims(icty,1:3,:,5) = pTB1;           % averted proportion
    Gestims(icty,1:3,:,6) = RdNumAbiUse;    % averted number
    Gestims(icty,:,:,7)   = AlTBInc2020;
    Gestims(icty,:,:,8)   = AlTBInc2035;
    Gestims(icty,:,:,9)   = DrTBInc2020;
    Gestims(icty,:,:,10)  = DrTBInc2035;

    % get summary indactors for plotting
    Prcestims(icty,:,:,1) = NavtInc/1e3;    % unit: thousands
    Prcestims(icty,:,:,2) = NavtMor/1e3;    % unit: thousands
    Prcestims(icty,:,:,3) = PavtInc*100;    % unit: %
    Prcestims(icty,:,:,4) = PavtMor*100;    % unit: %
    Prcestims(icty,:,:,5) = NavtTx1/1e3;    % unit: thousands
    Prcestims(icty,:,:,6) = NavtTx2/1e3;    % unit: thousands
    Prcestims(icty,:,:,7) = RCE*100;        % unit: %
    Prcestims(icty,:,:,8) = RAB/1e3;        % unit: thousands
    Prcestims(icty,:,:,9) = PavtTx2*100;    % unit: %
    Prcestims(icty,:,:,10) = Avt_Vx2Tx*100; % unit: % (case, death, PM)

    clear Allsoln data NCumCase NCumDeath NCumPMTx1 NCumPMTx2 pTB1 RdNumAbiUse .....
        NavtInc NavtMor PavtInc PavtMor NavtTx2 PavtTx2 RCE RAB 

    fprintf('Summary data processed for %s\n', CtyISO3{icty})

end
