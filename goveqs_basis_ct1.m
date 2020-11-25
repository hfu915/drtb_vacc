%% goveqs_basis_ct1.m
% run main ODEs and demographical changes

function out = goveqs_basis_ct1(t, in, M, i, s, r, p)

invec = in(1:i.nstates);

% Normalise by populations                     
lamCh = M.lambda(1:2,:)*invec/sum(invec(s.ch));      
lamAd = M.lambda(3:4,:)*invec/sum(invec(s.ad));       
lamFinal = p.contact*[lamCh lamAd]';

% Assign time-varying betas for both DS- and DR-TB
r_betaRt = r.betaRt_fun(t);  % r.betaRt_fun = @(t)betaRt_prm^(t-p.betaRt_yr);
BetaLam = r_betaRt.*lamFinal;

allmat = M.lin + BetaLam(1,1)*M.nlin.DS.ch + BetaLam(1,2)*M.nlin.DR.ch ...
               + BetaLam(2,1)*M.nlin.DS.ad + BetaLam(2,2)*M.nlin.DR.ad;
out = allmat*invec;

% Implement deaths
morts = sum(M.mort,1)'.*invec;                         
out = out - morts;

% Implement births
births = (p.equpop == 1)*sum(sum(morts)) + (p.equpop == 0)*r.birth*sum(invec);
% births = r.birth*sum(invec);
out(i.U.ch.v0) = out(i.U.ch.v0)+births*(1-p.vacc_birthday(1));
out(i.U.ch.v1) = out(i.U.ch.v1)+births*p.vacc_birthday(1);

% Get the auxiliaries
out(i.aux.inc)  = M.incdout*invec;
out(i.aux.ntfsec) = M.ntfsec*invec;
out(i.aux.newinf) = i.agg.newinf*(i.sel.newinf.*allmat)*invec;
out(i.aux.vdose) = sum((i.sel.vdose.*allmat)*invec);

tmp = M.mort([2,3],:)*invec;
out(i.aux.mort) = [sum(tmp) tmp'];  % All-TB, DS-TB, DR-TB 