function [Fit1Lorentzian,resnorm,exitflag,output,Lorentz1_1,cs_WB_1L] = fit1Lorentzian(x_data, y_data, initParams, x_fit, pH)

% fit the data to lorenzians
F_1_lorentzians = @(x,xdata)x(1) .* ( 1/pi .* x(2)./2 ./ ( (xdata-x(3)).^2 +(x(2)./2).^2 ) ) ;
% x(1) = amp
% x(2) = sigma [ppm]
% x(3) = x0 (center) [ppm]

[fit1Params,resnorm,resid,exitflag,output,lambda,Jacobian] = lsqcurvefit(F_1_lorentzians,initParams,x_data,y_data);

ci = nlparci(fit1Params,resid,'jacobian',Jacobian);

Fit1Lorentzian = F_1_lorentzians(fit1Params,x_fit);

Lorentz1_1 = fit1Params(1) .* (1/pi .* fit1Params(2)/2 ./ ( (x_fit-fit1Params(3)).^2 +(fit1Params(2)/2).^2 ));

cs_WB_1L        = fit1Params(3);
cs_WB_2L_err    = (ci(3,2)-ci(3,1))/2;
cs_WB_1L        = [cs_WB_1L,cs_WB_2L_err];

end


