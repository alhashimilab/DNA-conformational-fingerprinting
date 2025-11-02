function [Fit3Lorentzians,resnorm,exitflag,output,Lorentz1_3,Lorentz2_3,Lorentz3_3,p_WB_3L,p_Anion_3L,p_ES_3L,pKa_3L,cs_WB_3L,cs_Anion_3L,dcs_Anion_3L,cs_ES_3L] = fit3Lorentzians(x_data, y_data, initParams, x_fit, pH)

% fit the data to lorenzians
F_3_lorentzians = @(x,xdata)x(1) .* ( 1/pi .* x(2)./2 ./ ( (xdata-x(3)).^2 +(x(2)./2).^2 ) ) + x(4) .* ( 1/pi .* x(5)./2 ./ ( (xdata-x(6)).^2 +(x(5)./2).^2 ) ) + x(7) .* ( 1/pi .* x(8)./2 ./ ( (xdata-x(9)).^2 +(x(8)./2).^2 ) );
% x(1), x(4), x(7) = amp
% x(2), x(5), x(8) = sigma [ppm]
% x(3), x(6), x(9) = x0 (center) [ppm]

[fit3Params,resnorm,resid,exitflag,output,lambda,Jacobian] = lsqcurvefit(F_3_lorentzians,initParams,x_data,y_data);

ci = nlparci(fit3Params,resid,'jacobian',Jacobian);

Fit3Lorentzians = F_3_lorentzians(fit3Params,x_fit);

Lorentz1_3 = fit3Params(1) .* (1/pi .* fit3Params(2)/2 ./ ( (x_fit-fit3Params(3)).^2 +(fit3Params(2)/2).^2 ));
Lorentz2_3 = fit3Params(4) .* (1/pi .* fit3Params(5)/2 ./ ( (x_fit-fit3Params(6)).^2 +(fit3Params(5)/2).^2 ));
Lorentz3_3 = fit3Params(7) .* (1/pi .* fit3Params(8)/2 ./ ( (x_fit-fit3Params(9)).^2 +(fit3Params(8)/2).^2 ));

% calculate populations and pKa
p_WB_amp_3    = (fit3Params(4)/(fit3Params(1)+fit3Params(4)+fit3Params(7)))*100;
p_Anion_amp_3 = (fit3Params(1)/(fit3Params(1)+fit3Params(4)+fit3Params(7)))*100;
p_ES_amp_3    = (fit3Params(7)/(fit3Params(1)+fit3Params(4)+fit3Params(7)))*100;

amp_WB_3_err    = (ci(4,2)-ci(4,1))/2;
amp_Anion_3_err = (ci(1,2)-ci(1,1))/2;
amp_ES_3_err    = (ci(7,2)-ci(7,1))/2;

p_WB_amp_3_err = p_WB_amp_3*sqrt((amp_WB_3_err/fit3Params(4))^2+...
((amp_Anion_3_err)^2+(amp_WB_3_err)^2+(amp_ES_3_err)^2)/((fit3Params(1)+fit3Params(4)+fit3Params(7))^2));
p_Anion_amp_3_err = p_Anion_amp_3*sqrt((amp_Anion_3_err/fit3Params(1))^2+...
((amp_Anion_3_err)^2+(amp_WB_3_err)^2+(amp_ES_3_err)^2)/((fit3Params(1)+fit3Params(4)+fit3Params(7))^2));
p_ES_amp_3_err = p_ES_amp_3*sqrt((amp_ES_3_err/fit3Params(7))^2+...
((amp_Anion_3_err)^2+(amp_WB_3_err)^2+(amp_ES_3_err)^2)/((fit3Params(1)+fit3Params(4)+fit3Params(7))^2));

trapz1_3 = trapz(x_fit,Lorentz1_3);
trapz2_3 = trapz(x_fit,Lorentz2_3);
trapz3_3 = trapz(x_fit,Lorentz3_3);

p_WB_trapz_3    = (trapz2_3/(trapz1_3+trapz2_3+trapz3_3))*100;
p_Anion_trapz_3 = (trapz1_3/(trapz1_3+trapz2_3+trapz3_3))*100;
p_ES_trapz_3    = (trapz3_3/(trapz1_3+trapz2_3+trapz3_3))*100;

if p_WB_amp_3_err > std([p_WB_amp_3,p_WB_trapz_3])
    p_WB_3_err = p_WB_amp_3_err;
else
    p_WB_3_err = std([p_WB_amp_3,p_WB_trapz_3]);
end
    
if p_Anion_amp_3_err > std([p_Anion_amp_3,p_Anion_trapz_3])
    p_Anion_3_err = p_Anion_amp_3_err;
else
    p_Anion_3_err = std([p_Anion_amp_3,p_Anion_trapz_3]);
end

if p_ES_amp_3_err > std([p_ES_amp_3,p_ES_trapz_3])
    p_ES_3_err = p_ES_amp_3_err;
else
    p_ES_3_err = std([p_ES_amp_3,p_ES_trapz_3]);
end

p_WB_3L    = [mean([p_WB_amp_3,p_WB_trapz_3]),p_WB_3_err];
p_Anion_3L = [mean([p_Anion_amp_3,p_Anion_trapz_3]),p_Anion_3_err];
p_ES_3L    = [mean([p_ES_amp_3,p_ES_trapz_3]),p_ES_3_err];

pKa_3L(1,1)  = pH - log10(p_Anion_3L(1,1)./(100-p_Anion_3L(1,1)));
%pKa_3L(1,2)  = 1/(log(10))*sqrt((p_Anion_3L(1,2)).^2+(p_WB_3L(1,2)).^2+(p_ES_3L(1,2)).^2);
pKa_3L(1,2) = 1/(log(10))*sqrt((p_Anion_3L(1,2)./p_Anion_3L(1,1)).^2+((p_Anion_3L(1,2))./(100-p_Anion_3L(1,1))).^2);

cs_WB_3L        = fit3Params(6);
cs_WB_3L_err    = (ci(6,2)-ci(6,1))/2;
cs_WB_3L        = [cs_WB_3L,cs_WB_3L_err];

cs_Anion_3L     = fit3Params(3);
cs_Anion_3L_err = (ci(3,2)-ci(3,1))/2;
cs_Anion_3L     = [cs_Anion_3L,cs_Anion_3L_err];

dcs_Anion_3L    = cs_Anion_3L - cs_WB_3L;

cs_ES_3L        = fit3Params(9);
cs_ES_3L_err    = (ci(9,2)-ci(9,1))/2;
cs_ES_3L        = [cs_ES_3L,cs_ES_3L_err];
