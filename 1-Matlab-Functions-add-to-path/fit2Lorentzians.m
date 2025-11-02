function [Fit2Lorentzians,resnorm,exitflag,output,Lorentz1_2,Lorentz2_2,p_WB_2L,p_Anion_2L,pKa_2L,cs_WB_2L,cs_Anion_2L,dcs_Anion_2L] = fit2Lorentzians(x_data, y_data, initParams, x_fit, pH)

% fit the data to lorenzians
F_2_lorentzians = @(x,xdata)x(1) .* ( 1/pi .* x(2)./2 ./ ( (xdata-x(3)).^2 +(x(2)./2).^2 ) ) + x(4) .* ( 1/pi .* x(5)./2 ./ ( (xdata-x(6)).^2 +(x(5)./2).^2 ) );
% x(1), x(4) = amp
% x(2), x(5) = sigma [ppm]
% x(3), x(6) = x0 (center) [ppm]

[fit2Params,resnorm,resid,exitflag,output,lambda,Jacobian] = lsqcurvefit(F_2_lorentzians,initParams,x_data,y_data);

ci = nlparci(fit2Params,resid,'jacobian',Jacobian);

Fit2Lorentzians = F_2_lorentzians(fit2Params,x_fit);

Lorentz1_2 = fit2Params(1) .* (1/pi .* fit2Params(2)/2 ./ ( (x_fit-fit2Params(3)).^2 +(fit2Params(2)/2).^2 ));
Lorentz2_2 = fit2Params(4) .* (1/pi .* fit2Params(5)/2 ./ ( (x_fit-fit2Params(6)).^2 +(fit2Params(5)/2).^2 ));


% calculate populations and pKa
p_WB_amp_2    = (fit2Params(4)/(fit2Params(1)+fit2Params(4)))*100;
p_Anion_amp_2 = (fit2Params(1)/(fit2Params(1)+fit2Params(4)))*100;

amp_WB_2_err    = (ci(4,2)-ci(4,1))/2;
amp_Anion_2_err = (ci(1,2)-ci(1,1))/2;

p_WB_amp_2_err = p_WB_amp_2*sqrt((amp_WB_2_err/fit2Params(4))^2+...
((amp_Anion_2_err)^2+(amp_WB_2_err)^2)/((fit2Params(1)+fit2Params(4))^2));
p_Anion_amp_2_err = p_Anion_amp_2*sqrt((amp_Anion_2_err/fit2Params(1))^2+...
((amp_Anion_2_err)^2+(amp_WB_2_err)^2)/((fit2Params(1)+fit2Params(4))^2));

trapz1_2 = trapz(x_fit,Lorentz1_2);
trapz2_2 = trapz(x_fit,Lorentz2_2);

p_WB_trapz_2    = (trapz2_2/(trapz1_2+trapz2_2))*100;
p_Anion_trapz_2 = (trapz1_2/(trapz1_2+trapz2_2))*100;


if p_WB_amp_2_err > std([p_WB_amp_2,p_WB_trapz_2])
    p_WB_2_err = p_WB_amp_2_err;
else
    p_WB_2_err = std([p_WB_amp_2,p_WB_trapz_2]);
end
    
if p_Anion_amp_2_err > std([p_Anion_amp_2,p_Anion_trapz_2])
    p_Anion_2_err = p_Anion_amp_2_err;
else
    p_Anion_2_err = std([p_Anion_amp_2,p_Anion_trapz_2]);
end

p_WB_2L    = [mean([p_WB_amp_2,p_WB_trapz_2]),p_WB_2_err];
p_Anion_2L = [mean([p_Anion_amp_2,p_Anion_trapz_2]),p_Anion_2_err];

pKa_2L(1,1)  = pH - log10(p_Anion_2L(1,1)./(100-p_Anion_2L(1,1)));
%pKa_2L(1,2)  = 1/(log(10))*sqrt((p_Anion_2L(1,2)).^2+(p_WB_2L(1,2)).^2);
pKa_2L(1,2) = 1/(log(10))*sqrt((p_Anion_2L(1,2)./p_Anion_2L(1,1)).^2+((p_Anion_2L(1,2))./(100-p_Anion_2L(1,1))).^2);

cs_WB_2L        = fit2Params(6);
cs_WB_2L_err    = (ci(6,2)-ci(6,1))/2;
cs_WB_2L        = [cs_WB_2L,cs_WB_2L_err];

cs_Anion_2L     = fit2Params(3);
cs_Anion_2L_err = (ci(3,2)-ci(3,1))/2;
cs_Anion_2L     = [cs_Anion_2L,cs_Anion_2L_err];

dcs_Anion_2L = cs_Anion_2L - cs_WB_2L;


