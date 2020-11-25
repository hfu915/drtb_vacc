%% make_model_ct1.m
% set up matrices for transitions between compartments

function M = make_model_ct1(p, r, i, s, gps, vacc_inds)

nostates = i.nstates;
m = zeros(nostates);
minc = zeros(nostates); 
mntfsec = zeros(2,nostates);
macq = zeros(nostates);

for ia = 1:length(gps.ages)
    age = gps.ages{ia};
    
    for iv = 1:length(gps.vaccs)
        vacc = gps.vaccs{iv};
        
        for istr = 1:length(gps.strains)
            strain = gps.strains{istr};
            
            getadd = @(st) i.(st).(age).(vacc).(strain);
            
            % --- Get the linear rates ------------------------------------------------
            L0    = getadd('L0');
            L1    = getadd('L1');
            I     = getadd('I');
            Dxs   = getadd('Dx');
            Txs   = getadd('Tx');
            Tx2s  = getadd('Tx2');
            E     = getadd('E');
            Rlo   = getadd('Rlo');
            Rhi   = getadd('Rhi');
            R     = getadd('R');
            
            % --- Breakdown and reactivation
            sources = [L0, L1]; destin  = I;
            rates   = [r.breakdown(ia,iv), r.reactivation(ia,iv)];
            minc(destin, sources) = minc(destin, sources) + rates;
            
            % --- Stabilisation
            source = L0; destin = L1; rate = r.stabilisation;
            m(destin, source) = m(destin, source) + rate;
            
            % --- Primary careseeking, including access to public sector care
            source = I; destin = Dxs.pu; rate = r.careseeking*p.pu; 
            m(destin, source) = m(destin, source) + rate;
            
            source = I; destin = Dxs.pr; rate = r.careseeking*(1-p.pu);
            m(destin, source) = m(destin, source) + rate;
            
            for ip = 1:length(gps.sectors)
                prov = gps.sectors{ip};
                Dx = Dxs.(prov); Tx = Txs.(prov); Tx2 = Tx2s.(prov);
                
                % --- Diagnosis
                ismdr   = strcmp(strain, 'DR');
                pFLinit = p.Tx_initFL(ip)*(1 - ismdr*p.DR_rec(ip));
                pSLinit = p.Tx_initSL(ip)*ismdr*p.DR_rec(ip);
                p_ltfu  = 1-(pFLinit+pSLinit);
                
                source  = Dx;
                destins =          [Tx,       Tx2,      E];
                rates   = r.Dx*[pFLinit,  pSLinit,  p_ltfu];
                m(destins, source) = m(destins, source) + rates';
                
                mntfsec(ip,Dx) = mntfsec(ip,Dx) + r.Dx*(pFLinit+pSLinit);
                
                % --- FL Treatment
                pFLcure = p.cureFL(ip)*(1-ismdr);
                pSLtran = p.SL_trans(ip)*ismdr;
                               
                source  = Tx;
                destins = [Rlo              , Tx2                          , E                                , Rhi                            ];
                rates   = [r.txFLcom*pFLcure, r.txFLcom*(1-pFLcure)*pSLtran, r.txFLcom*(1-pFLcure)*(1-pSLtran), r.self_cure*ismdr+r.txFLdef(ip)];
                m(destins, source) = m(destins, source) + rates';
                
				% --- DR acquisition from FL treatment
                source = Tx;
                destin = i.Tx.(age).(vacc).DR.(prov);
                rate   = r.DR_acqu*(1-ismdr);
                macq(destin, source) = macq(destin, source) + rate;
					
                % --- SL Treatment
                source  = Tx2;
                destins = [Rlo                   , E                                   ];
                rates   = [r.txSLcom*p.cureSL(ip), r.txSLcom*(1-p.cureSL(ip))+r.txSLdef];
                m(destins, source) = m(destins, source) + rates';
                
            end
            
            % --- Secondary careseeking
            source = E; destin = Dxs.pu; rate = r.careseeking2*p.pu;
            m(destin, source) = m(destin, source) + rate;
            
            source = E; destin = Dxs.pr; rate = r.careseeking2*(1-p.pu);
            m(destin, source) = m(destin, source) + rate;
            
            % --- Relapse
            sources = [Rlo Rhi R];
            destin  = I;
            rates   = r.relapse;
            minc(destin, sources) = minc(destin, sources) + rates;
            
            sources = [Rlo Rhi];
            destin  = R;
            rates   = 0.5;
            m(destin, sources) = m(destin, sources) + rates;
            
            % --- Self cure
            sources = [I Dx E]; 
            destin  = Rhi;
            rates   = r.self_cure;
            m(destin, sources) = m(destin, sources) + rates;
            
        end
    end
end


% % --- Age only, without vaccination
% inds = sub2ind([nostates,nostates], s.ad, s.ch);
% m(inds) = m(inds) + r.ageing;

% --- Ageing and vaccination
getinds_vac = @(age,vacc) intersect(intersect(s.(age),s.(vacc)),vacc_inds);

% --- Catch-up campaign
inds = sub2ind([nostates,nostates], getinds_vac('ch','v1'), getinds_vac('ch','v0'));
m(inds) = m(inds) + r.vacc_catchup(1);

inds = sub2ind([nostates,nostates], getinds_vac('ad','v1'), getinds_vac('ad','v0'));
m(inds) = m(inds) + r.vacc_catchup(2);
 
% --- Ageing and regular vaccination on 15th birthday
inds = sub2ind([nostates,nostates], getinds_vac('ad','v1'), getinds_vac('ch','v0'));
m(inds) = m(inds) + r.ageing*p.vacc_birthday(2);

inds = sub2ind([nostates,nostates], getinds_vac('ad','v0'), getinds_vac('ch','v0'));
m(inds) = m(inds) + r.ageing*(1-p.vacc_birthday(2));

novacc_inds = setdiff(1:nostates,vacc_inds);
getinds_novac = @(age,vacc) intersect(intersect(s.(age),s.(vacc)),novacc_inds);

inds = sub2ind([nostates,nostates], getinds_novac('ad','v0'), getinds_novac('ch','v0'));
m(inds) = m(inds) + r.ageing;

inds = sub2ind([nostates,nostates], intersect(s.ad,s.v1), intersect(s.ch,s.v1));  % those with effective immunity
m(inds) = m(inds) + r.ageing;

inds = sub2ind([nostates,nostates], intersect(s.ad,s.v2), intersect(s.ch,s.v2));  % those with waned immunity    
m(inds) = m(inds) + r.ageing;

% --- Vaccine immunity waning
inds = sub2ind([nostates,nostates], s.v2, s.v1);
m(inds) = m(inds) + r.waning;

Mlin = m + minc + macq;
M.lin = sparse(Mlin - diag(sum(Mlin,1)));

% --- Get notifications in the sectors
M.ntfsec = sparse(mntfsec);

% --- Get the nonlinear rates ---------------------------------------------

% --- Allocating transmission
for istr = 1:length(gps.strains)
    m = zeros(nostates);
    for ia = 1:length(gps.ages)
        age = gps.ages{ia};
        for iv = 1:length(gps.vaccs)
            vacc = gps.vaccs{iv};
            
            strain = gps.strains{istr};
            U = i.U.(age).(vacc); L0 = i.L0.(age).(vacc).(strain); 
            imms = intersect(intersect([s.L0, s.L1, s.Rlo, s.Rhi, s.R], s.(age)),s.(vacc));
            
            m(L0, [U, imms]) = 1;
            
        end
        m(:,[s.L0, s.L1, s.Rlo, s.Rhi, s.R]) = m(:,[s.L0, s.L1, s.Rlo, s.Rhi, s.R])*p.red_sus;
        M.nlin.(strain).(age) = sparse(m - diag(sum(m,1)));
    end
end


% --- Getting force-of-infection
m = zeros(4,nostates);
m(1,intersect(intersect(s.infectious,s.DS),s.ch)) = r.beta(1)*p.red_infectCh;
m(2,intersect(intersect(s.infectious,s.DR),s.ch)) = r.beta(2)*p.red_infectCh;
m(3,intersect(intersect(s.infectious,s.DS),s.ad)) = r.beta(1);
m(4,intersect(intersect(s.infectious,s.DR),s.ad)) = r.beta(2);
M.lambda = sparse(m);

% --- Get incidence rates: all-form, DS(ch), DS(ad), DR(ch), DR(ad), DR acq
m = zeros(5,nostates);
m(1,[s.I s.Tx]) = 1;    
m(2,intersect(intersect(s.ch,s.DS),s.I)) = 1;
m(3,intersect(intersect(s.ad,s.DS),s.I)) = 1;
m(4,intersect(intersect(s.ch,s.DR),[s.I s.Tx])) = 1;
m(5,intersect(intersect(s.ad,s.DR),[s.I s.Tx])) = 1;
M.incdout = sparse(m*(minc+macq));
M.incdout(6,:) = sum(macq,1);

% --- Get the mortality rates: all-cause, DS, DR
m = zeros(3,nostates);
m(1,setdiff(s.ch,[s.Tx s.Tx2]))            = r.mort(1);   
m(1,setdiff(s.ad,[s.Tx s.Tx2]))            = r.mort(2);     
m(2,intersect([s.I s.Dx s.E],s.DS))        = r.mort_TB;
m(2,intersect(intersect(s.Tx,s.pu),s.DS))  = r.txFLdie(1);  
m(2,intersect(intersect(s.Tx,s.pr),s.DS))  = r.txFLdie(2);                        
m(3,intersect([s.I s.Dx s.E],s.DR))        = r.mort_TB;
m(3,intersect(intersect(s.Tx,s.pu),s.DR))  = r.txFLdie(1);   
m(3,intersect(intersect(s.Tx,s.pr),s.DR))  = r.txFLdie(2);                          
m(3,intersect(s.Tx2,s.DR))                 = r.txSLdie;                       
M.mort = sparse(m);     
