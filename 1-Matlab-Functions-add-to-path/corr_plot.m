function [r,Rsq,Rsq_det,rms,a_rms_per,bestfit,x2,inBetween,slope,intercept,slope_err,intercept_err,p] = corr_plot(x_data,y_data,x_line,conf_int_per)

[r_corr,p_corr] = corrcoef(x_data,y_data);
r          = r_corr(1,2);
p          = p_corr(1,2);

linreg     = fit(x_data,y_data,'poly1');
mdl        = fitlm(x_data,y_data);
Rsq        = mdl.Rsquared.Ordinary;
slope      = mdl.Coefficients.Estimate(2);
intercept  = mdl.Coefficients.Estimate(1);
slope_err  = mdl.Coefficients.SE(2);
intercept_err = mdl.Coefficients.SE(1);
rms        = sqrt(1/(length(x_data)).*(sum((y_data - x_data).^2)));

a_rms_per  = 100*rms/(max(y_data) - min(y_data));

y_pred     = x_data;
y_actual   = y_data;
SSR        = sum((y_pred - y_actual).^2);
TSS        = sum(((y_actual - mean(y_actual)).^2));
Rsq_det    = 1 - SSR/TSS;

ci         = confint(linreg,conf_int_per);

ci_line_1  = ci(1,1).*x_line + ci(1,2); % p1_lb*x + p2_lb
ci_line_2  = ci(2,1).*x_line + ci(2,2); % p1_ub*x + p2_ub

ci_line_3  = ci(2,1).*x_line + ci(1,2); % p1_ub*x + p2_lb
ci_line_4  = ci(1,1).*x_line + ci(2,2); % p1_lb*x + p2_ub

ci_max     = max([ci_line_1;ci_line_2;ci_line_3;ci_line_4]);
ci_min     = min([ci_line_1;ci_line_2;ci_line_3;ci_line_4]);

bestfit    = linreg(x_line);

x2         = [x_line, fliplr(x_line)];
inBetween  = [ci_min, fliplr(ci_max)];

end