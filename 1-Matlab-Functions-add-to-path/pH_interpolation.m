function [params_pH_tbl] = pH_interpolation(seq, pH_for_calc, pH_for_calc_errs, ...
    pHs, pKas, pKas_err, ...
    pGS, pGS_err, pT, pT_err, pA, pA_err, ...
    kGStoES1, kGStoES1_err, kES1toGS, kES1toGS_err, ...
    kGStoES2, kGStoES2_err, kES2toGS, kES2toGS_err, ...
    kES2toES1, kES2toES1_err, kES1toES2, kES1toES2_err, ...
    slopes_GStoES2, slopes_minor)

%make calculations with fractions (not percents)

pT_pH      = pT./100;
pT_pH_err  = pT_err./100;

pA_pH      = (10.^(pH_for_calc - pKas))./(1+10.^(pH_for_calc - pKas));
pA_pH_err  = pA_pH.*sqrt((pKas_err.*log(10)).^2+((10.^(pH_for_calc - pKas).*(pKas_err.*log(10)))./(1+10.^(pH_for_calc - pKas))).^2);

pA_pH(pHs==pH_for_calc,:)     = pA(pHs==pH_for_calc,:)./100;
pA_pH_err(pHs==pH_for_calc,:) = pA_err(pHs==pH_for_calc,:)./100;

pGS_pH     = 1 - pA_pH - pT_pH;
pGS_pH_err = sqrt((pA_pH_err).^2+(pT_pH_err).^2);

pGS_pH(pHs==pH_for_calc,:)     = pGS(pHs==pH_for_calc,:)./100;
pGS_pH_err(pHs==pH_for_calc,:) = pGS_err(pHs==pH_for_calc,:)./100;

kGStoES1_pH               = kGStoES1;
kGStoES1_pH_err           = kGStoES1_err;
kES1toGS_pH               = kES1toGS;
kES1toGS_pH_err           = kES1toGS_err;


kGStoES2_for_pH_ln        = log(kGStoES2);
kGStoES2_for_pH_ln_err    = kGStoES2_err./kGStoES2;
inter_GStoES2_for_pH      = kGStoES2_for_pH_ln - slopes_GStoES2(1,3).*pHs;
inter_GStoES2_for_pH_err  = sqrt((kGStoES2_for_pH_ln_err).^2+(pHs).^2.*(slopes_GStoES2(2,3)).^2);
kGStoES2_pH_ln            = slopes_GStoES2(1,3).*pH_for_calc + inter_GStoES2_for_pH;
kGStoES2_pH_ln_err        = sqrt((inter_GStoES2_for_pH_err).^2+(pH_for_calc).^2.*(slopes_GStoES2(2,3)).^2);
kGStoES2_pH               = exp(kGStoES2_pH_ln);
kGStoES2_pH_err           = kGStoES2_pH.*kGStoES2_pH_ln_err;

kGStoES2_pH(pHs==pH_for_calc,:)     = kGStoES2(pHs==pH_for_calc,:);
kGStoES2_pH_err(pHs==pH_for_calc,:) = kGStoES2_err(pHs==pH_for_calc,:);

kES2toGS_pH               = kES2toGS;
kES2toGS_pH_err           = kES2toGS_err;


kES2toES1_for_pH_ln       = log(kES2toES1);
kES2toES1_for_pH_ln_err   = kES2toES1_err./kES2toES1;
inter_ES2toES1_for_pH     = kES2toES1_for_pH_ln - slopes_minor(1,1).*pHs;
inter_ES2toES1_for_pH_err = sqrt((kES2toES1_for_pH_ln_err).^2+(pHs).^2.*(slopes_minor(2,1)).^2);
kES2toES1_pH_ln           = slopes_minor(1,1).*pH_for_calc + inter_ES2toES1_for_pH;
kES2toES1_pH_ln_err       = sqrt((inter_ES2toES1_for_pH_err).^2+(pH_for_calc).^2.*(slopes_minor(2,1)).^2);
kES2toES1_pH              = exp(kES2toES1_pH_ln);
kES2toES1_pH_err          = kES2toES1_pH.*kES2toES1_pH_ln_err;

kES2toES1_pH(pHs==pH_for_calc,:)     = kES2toES1(pHs==pH_for_calc,:);
kES2toES1_pH_err(pHs==pH_for_calc,:) = kES2toES1_err(pHs==pH_for_calc,:);

kES1toES2_for_pH_ln       = log(kES1toES2);
kES1toES2_for_pH_ln_err   = kES1toES2_err./kES1toES2;
inter_ES1toES2_for_pH     = kES1toES2_for_pH_ln - slopes_minor(1,2).*pHs;
inter_ES1toES2_for_pH_err = sqrt((kES1toES2_for_pH_ln_err).^2+(pHs).^2.*(slopes_minor(2,2)).^2);
kES1toES2_pH_ln           = slopes_minor(1,2).*pH_for_calc + inter_ES1toES2_for_pH;
kES1toES2_pH_ln_err       = sqrt((inter_ES1toES2_for_pH_err).^2+(pH_for_calc).^2.*(slopes_minor(2,2)).^2);
kES1toES2_pH              = exp(kES1toES2_pH_ln);
kES1toES2_pH_err          = kES1toES2_pH.*kES1toES2_pH_ln_err;

kES1toES2_pH(pHs==pH_for_calc,:)     = kES1toES2(pHs==pH_for_calc,:);
kES1toES2_pH_err(pHs==pH_for_calc,:) = kES1toES2_err(pHs==pH_for_calc,:);


kES1_pH     = kGStoES1_pH + kES1toGS_pH;
kES1_pH_err = sqrt(kGStoES1_pH_err.^2 + kES1toGS_pH_err.^2);

kES2_pH     = kGStoES2_pH + kES2toGS_pH;
kES2_pH_err = sqrt(kGStoES2_pH_err.^2 + kES2toGS_pH_err.^2);

kMinor_pH     = kES1toES2_pH + kES2toES1_pH;
kMinor_pH_err = sqrt(kES1toES2_pH_err.^2 + kES2toES1_pH_err.^2);

%return parameters (populations are fractions, not percents)
params_pH = [pGS_pH,pGS_pH_err,pT_pH,pT_pH_err,pA_pH,pA_pH_err,...
    kES1_pH,kES1_pH_err,kES2_pH,kES2_pH_err,kMinor_pH,kMinor_pH_err,...
    kGStoES1_pH,kGStoES1_pH_err,kES1toGS_pH,kES1toGS_pH_err,...
    kGStoES2_pH,kGStoES2_pH_err,kES2toGS_pH,kES2toGS_pH_err,...
    kES1toES2_pH,kES1toES2_pH_err,kES2toES1_pH,kES2toES1_pH_err];

params_titles = {'pGS_pH','pGS_pH_err','pT_pH','pT_pH_err','pA_pH','pA_pH_err',...
    'kES1_pH','kES1_pH_err','kES2_pH','kES2_pH_err','kMinor_pH','kMinor_pH_err',...
    'kGStoES1_pH','kGStoES1_pH_err','kES1toGS_pH','kES1toGS_pH_err',...
    'kGStoES2_pH','kGStoES2_pH_err','kES2toGS_pH','kES2toGS_pH_err',...
    'kES1toES2_pH','kES1toES2_pH_err','kES2toES1_pH','kES2toES1_pH_err'};

%params_pH_to_save = [params_titles;num2cell(params_pH)];
orig_pH = pHs;

params_pH_tbl = array2table(params_pH);
params_pH_tbl.Properties.VariableNames = params_titles;
params_pH_tbl = addvars(params_pH_tbl,seq,'Before','pGS_pH');
params_pH_tbl = addvars(params_pH_tbl,pH_for_calc,'Before','pGS_pH');
params_pH_tbl = addvars(params_pH_tbl,pH_for_calc_errs,'Before','pGS_pH');
params_pH_tbl = addvars(params_pH_tbl,orig_pH,'Before','pGS_pH');


