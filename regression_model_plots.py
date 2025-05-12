##########################################################
## Compute current and pre-industrial conditions, conceptually
## add alkalinity (in the form of NaOH and Na2CO3) 
##########################################################
## Author: H. van de Mortel
## Version: 1.0
## Maintainer: H. van de Mortel
## Email: vandemortel.hanna@gmail.com
## Using funtions by Greg Pelletier
##########################################################

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
from sklearn.preprocessing import PolynomialFeatures
from scipy.stats import zscore
from numpy import exp, linspace
import lmfit as fit
from lmfit.models import ExpressionModel
from scipy import stats
from statsmodels.tools.eval_measures import rmse
import re
import matplotlib.colors as mcolors
import matplotlib.ticker as ticker

#%%
filepath = "C:/Users/hanna/Documents/GitHub/OAE_mitigation_of_OA/"

#Import compiled OA bio data
df1 = pd.read_excel(filepath + "master_data.xlsx", 
                    sheet_name='Sheet1')
#Import added alkalinity data
df2 = pd.read_excel(filepath + "output_TA_addition.xlsx", 
                    sheet_name='Sheet1')

#Define variables for simplicity
TA = 'AT [µmol/kg] (compiled)'
DIC = 'DIC [µmol/kg] (compiled)'
species2 = 'Species (compiled)'

rate1 = 'Calc rate [mmol/g/hr] (compiled)'
rate2 = 'Calc rate [mmol/m**2/hr] (compiled)'
rate3 = 'Calc rate [mmol/m**3/hr] (compiled)'
rate4 = 'Calc rate [mmol/#/hr] (compiled)'
rate5 = 'Calc rate [mmol/hr] (compiled)'
rate6 = 'Calc rate [mmol/cm**2] (compiled)'
rate7 = 'Calc rate [%/hr] (compiled)'

#Convert numeric columns to numeric type and handle non-numeric values
df1[DIC] = pd.to_numeric(df1[DIC], errors='coerce')
df1[TA] = pd.to_numeric(df1[TA], errors='coerce')

#Calculate TA-DIC column
nan_mask = np.logical_or(np.isnan(df1[DIC]), np.isnan(df1[TA]))
df1['TA-DIC'] = np.where(nan_mask, np.nan, df1[TA] - df1[DIC])

#Group by species for plotting
grouped = df1.groupby(species2)

#Define a colormap and the number of unique colors for plotting from different studies
num_unique_colors = 3
orange_colors = ['#8da0cb', '#66c2a5', 'xkcd:blush']  # Shades of orange
cmap = mcolors.LinearSegmentedColormap.from_list("custom_orange", orange_colors, N=3)
orange_colors = [cmap(i / (3 - 1)) for i in range(3)]

#Choose colors for +TA data
cmap_NaOH =[
    "#b2e0d8",  # Light teal
    "#66a69f",  # Medium-light teal
    "#2d7a7b",  # Deeper teal
    "#0d403d"   # Darkest teal
]


#%% Make lists for best fits dataframe
best_fits_studies = [] #studies
best_fits_n = [] #number of datapoints
best_fits_groups = [] #functional group
best_fits_species = [] #species
best_fits_rate = [] #rate unit
best_fits = [] #response (linear, parabolic or exponential)
best_fits_slope =[] #slope for linear
best_fits_intercept = []
best_fits_exponent = []
best_fits_p = [] #p-value
best_fits_r2 = [] #R^2
best_fits_rmse = [] #RMSE

best_fits_pre_industrial= []
best_fits_baseline_pi = []

best_fits_baseline =[]
best_fits_current = []

best_fits_NaOH_50_uneq = []
best_fits_NaOH_50_eq = []

best_fits_Na2CO3_50_uneq = []
best_fits_Na2CO3_50_eq = []

best_fits_TA_pi =[]
best_fits_DIC_pi = []

best_fits_TA_baseline =[]
best_fits_DIC_baseline = []

best_fits_lower = []
best_fits_upper = []
best_fits_lower_pi = []
best_fits_upper_pi = []

best_fits_x_range = []

#Make lists for thresholds
NaOH_threshold_species = []
NaOH_threshold_groups = []
NaOH_threshold_rate = []
NaOH_threshold_amount = []
NaOH_thresholds = []

Na2CO3_threshold_species = []
Na2CO3_threshold_groups = []
Na2CO3_threshold_rate = []
Na2CO3_threshold_amount = []
Na2CO3_thresholds = []

TA_ranges = []
DIC_ranges = []

#%%

# function to use lmfit to find the best-fit parameters of the 3-parameter exponential regression for threshold+ or threshold- response
def fit_exp3(x,y, nper):
    # - - -
    # function to calculate the best-fit nonlinear regression parameters 
    # of the 3-parameter exponential regression for threshold+ or threshold- response
    # using lmfit (https://lmfit.github.io//lmfit-py/index.html)
    # Greg Pelletier (gjpelletier@gmail.com)
    # - - -

    # INPUT
    # x = vector of observed x
    # y = vector ov observed y
    #
    # OUTPUT
    # param = vector of best-fit parameters of the nonlinear regression function

    # find y stats to use for intial parameter estimates
    ymin = np.nanmin(y)
    yn = np.percentile(y,nper)
    ymed = np.median(y)
    ymax = np.nanmax(y,0)
    
    mod = ExpressionModel('b1 + b2 * exp(b3 * x)')
    pars = fit.Parameters()

    # - - -
    # This is where you set the initial values of the regression parameters b1, b2, b3 in pars.
    # The initial guess of parameters set to ypred = b1 = ymed (np.median(y)) usually is the best guess for intial b1
    # Initial values of b2=0 and b3=0 usually works best
    # In other words, the following initial values of b1, b2, b3 as follows usually works:
    # b1 = ymed
    # b2 = 0
    # b3 = 0
    # Sometimes this does not work. 
    # You can tell that it does not work when the best fit is a linear straight line through the observations
    # instead of a curving exponential threshold line.
    # In that case, try using either b1=y25th or b1=y75th instead of b1=med
    # Usually one of those for initial b1 will find a better fit with a curving exponential thrshold line when the b1=ymed does not work
    # In this example, b1=ymed does not work, but b1=y25th does work
    # Uncomment one of the following three pars['b1'] lines. Try ymed first. If that does not work then try y25th and y75th, and pick whichever works best
    # - - -
    #if yb=='ymed':
    #    pars['b1'] = fit.Parameter(name='b1', value=ymed, min=-np.inf, max=np.inf)
    #elif yb=='y25th':
    #    pars['b1'] = fit.Parameter(name='b1', value=y25th, min=-np.inf, max=np.inf)
    #else:
    #    pars['b1'] = fit.Parameter(name='b1', value=y75th, min=-np.inf, max=np.inf)
        
    pars['b1'] = fit.Parameter(name='b1', value=yn, min=-np.inf, max=np.inf)
    pars['b2'] = fit.Parameter(name='b2', value=0, min=-np.inf, max=np.inf)
    pars['b3'] = fit.Parameter(name='b3', value=0, min=-np.inf, max=np.inf)

    # extract parameter values from the best fit output and make user-defined function for the regression
    out = mod.fit(y, pars, x=x)
    d = out.params
    b1 = d['b1'].value
    b2 = d['b2'].value 
    b3 = d['b3'].value 
    param = [b1, b2, b3]
    
    return out, param


# function to calculate the confidence interval or prediction interval
def prediction_interval(f,param,xpred,x,y,alpha,prediction):
    # - - -
    # function to calculate the confidence interval or prediction interval
    # for any user-defined regression function.
    # Greg Pelletier (gjpelletier@gmail.com)
    # - - -

    # INPUT
    # f = user-defined regression @ function to predict y given inputs if param and x values (xpred and x)
    # 	For example, if using the 3-parameter nonlinear regression exponential threshold function, then 
    # 	f = lambda param,xval : param[0] + param[1] * exp(param[2] * xval)
    # param = vector of best-fit parameters of the regression function
    # xpred = vector of x values to evaluate predicted y values (e.g. xpred=linspace(min(x),max(x))
    # x = vector of observed x
    # y = vector ov observed y
    # alpha = probability value for the prediction interval (e.g. alpha=0.05 is the 95# prediction interval)
    # prediction = True or False, where True= prediction interval of additional measurements, and False= confidence iterval of original observations
    #
    # OUTPUT
    # ypred = predicted y values at xpred
    # ciLo = lower prediction interval or confidence interval for each value in xpred
    # ciUp = upper prediction interval or confidence interval for each value in xpred
    # p-value of the regression F-test comparing MSregression/MSresidual

    # use 2-tailed t statistic
    pLo = alpha/2
    pUp = 1-alpha/2
    nobs = np.size(x)
    nu = nobs-np.size(param)    # degrees of freedom = number of samples - number of parameters

    # critical t values for the upper and lower CI
    tcrit = stats.t.ppf([pLo, pUp], nu)
    
    # predicted ypred at xpred and yhat at x
    ypred = f(param,xpred)
    yhat = f(param,x)
    
    # residual sum of squares
    SSresidual = 0
    for i in range(nobs):
        SSresidual = SSresidual + (y[i] - yhat[i]) ** 2
    # residual mean square
    MSresidual = SSresidual / (np.size(x)-np.size(param))
    
    # syx = standard error of the estimate 
    syx = np.sqrt(MSresidual)

    # mean value of observed x
    xmean = np.mean(x)

    # sum of squares for x
    SSx = np.sum(x **2) - np.sum(x) **2 / nobs

    # calculate F-statistic and p-value of the regression
    SStotal = np.sum(y **2) - np.sum(y) **2 / nobs
    SSregression = SStotal - SSresidual
    # MSregression = SSregression
    MSregression = SSregression / (np.size(param)-1)
    Fstat = MSregression / MSresidual
    dfn = np.size(param) - 1                # df numerator = degrees of freedom for model = number of model parameters - 1
    dfd = np.size(x) - np.size(param)       # df denomenator = degrees of freedom of the residual = nobs - nparam
    p = 1-stats.f.cdf(Fstat, dfn, dfd)      # p-value of F test statistic 

    # calculate rsquared (ordinary and adjusted)
    rsquared = SSregression / SStotal                                           # ordinary rsquared
    adj_rsquared = 1-(1-rsquared)*(np.size(x)-1)/(np.size(x)-np.size(param)-1)  # adjusted rsquared

    # calculate sqrterm
    npred = np.size(xpred)
    sqrterm = np.empty(npred)
    for i in range(npred):
        if prediction:
            # prediction interval for additional observations
            sqrterm[i] = np.sqrt(1 + 1/nobs + (xpred[i] - xmean)**2 / SSx)
        else:
            # confidence interval of original observations
            sqrterm[i] = np.sqrt(1/nobs + (xpred[i] - xmean)**2 / SSx)

    ciLo = ypred + tcrit[0] * syx * sqrterm
    ciUp = ypred + tcrit[1] * syx * sqrterm

    return ypred, ciLo, ciUp, p, rsquared

# - - -
# function to calculate the adjusted r-squared
def adjusted_rsquared(rsquared,n,p):
    # calculate the adjusted rsquared from input of the following:
    # rsquared
    # n = number of samples
    # p = number of parameters
    adj_rsquared = 1-(1-rsquared)*(n-1)/(n-p-1)
    return adj_rsquared


#%%
for rate in [rate1, rate2, rate3, rate4, rate5, rate6, rate7]:
    df1[rate] = pd.to_numeric(df1[rate], errors='coerce')

s_species, s_rate, pqs, pls, s_iyval, pes = list(), list(), list(), list(), list(), list()
perc = [10,25,50,75,90]

for species, group in grouped:
    for rate in [rate1, rate2, rate3, rate4, rate5, rate6, rate7]:                    
        rate_group = group[group[rate].notna()]
        rate_group[rate] = pd.to_numeric(rate_group[rate], errors='coerce')
            
        # Determine the number of unique 'File Name' values in the current species group
        unique_files = rate_group['Study'].nunique()
        unique_file_names = ''  # Initialize an empty string to hold the combined study names

        # Generate a list of colors based on the number of unique 'File Name' values
        colors = [cmap(i) for i in np.linspace(0, 1, num_unique_colors)[:unique_files]]
        
        ######## TA:DIC vs calc rate ########        
        #Create mask to identify NaN
        x_ = rate_group['TA-DIC']
        y_ = rate_group[rate]
        nan_mask = np.logical_or(pd.isnull(x_), pd.isnull(y_))
        x = x_[~nan_mask]
        y = y_[~nan_mask]
            
        # Check if the data contains valid values to perform regressions
        if len(y) > 3:
            #Remove outliers (up to 6 std deviations)
            z_scores = zscore(y)
            z_threshold = 6
            outlier_mask = np.abs(z_scores) < z_threshold
            x = x[outlier_mask]
            y = y[outlier_mask]
                    
            #Filter df2 & df3 based on species
            df2_ = df2[(df2['Species'] == species) & (df2['Rate'] == rate)]
            
            #Compute calc rate for new TA values
            # Define the column prefixes and suffixes
            prefixes = ['NaOH TA-DIC +']
            suffix = ' umol/kg'
            
            # Create the array for NaOH TA:DIC
            naoh_columns = [df2_[f'{prefixes[0]}{i}{suffix}'] for i in list(range(50, 201, 50))]
            
            # Concatenate both arrays into x1
            x1 = np.array(naoh_columns)
          
            labels =list(range(50, 201, 50))

            #Append study names
            for i, study in enumerate(rate_group['Study'][~nan_mask].unique()):
                if i > 0:
                     unique_file_names += ', '  
                unique_file_names += study
                
            #print(sum(x1)) # all > 0
            x_array = x.to_numpy()
            xpred_l = linspace(x.min(), x.max(), 100)
            
            # linear
            polynomial_features = PolynomialFeatures(degree=1)
            xp = polynomial_features.fit_transform(x_array.reshape(-1, 1))
            model_lin = sm.OLS(y, xp).fit()
            lp=model_lin.f_pvalue
            lrsqu=model_lin.rsquared
            lparams = model_lin.params
            fl = lambda lparams, xval : lparams[0] + lparams[1]*xval 
            ypred_lin = fl(lparams, x)
            lrmse = rmse(y, ypred_lin)
            
            # exponential
            ps, rsqs = list(), list()
            xmean = np.mean(x)

            for j in perc:
                out, param = fit_exp3(x,y, j)
                f = lambda param,xval : param[0] + param[1] * exp(param[2] * xval)

                # use function to find ypred, ciLo, ciUp, and p-value
                # on the range of the data
                ypred_l, ciLo, ciUp, p, rsq = prediction_interval(f,param,xpred_l,x.values,y.values,.1,True)
                ps.append(p)
                rsqs.append(rsq)
            xp = np.min(ps)
            mi = ps.index(xp)
            xqsqu = rsqs[mi]
            
            # recalculate exponential
            out, xparam = fit_exp3(x,y, perc[mi])
            fe = lambda xparam,xval : xparam[0] + xparam[1] * exp(xparam[2] * xval)
            # use function to find ypred, ciLo, ciUp, and p-value
            # on the range of the data
            ypred_exp = fe(xparam, x)
            ypred_exp_l, ciLo, ciUp, xp, rsq = prediction_interval(fe,xparam,xpred_l,x.values,y.values,.1,True)
            xrmse = rmse(y, ypred_exp)
                
            fig = plt.figure(figsize=(3,3),dpi=400)

            p_thr = 0.05
            if (lp<=p_thr) | (xp<=p_thr):                
                # plot data
                studies_handles = []
                study_labels = ['Experimental data from study A', 'Experimental data from study B', 'Experimental data from study C', 'Experimental data from study D']

                #Plot the data points with different colors for each 'File Name' value
                for i, study in enumerate(rate_group['Study'][~nan_mask].unique()):
                    mask = rate_group['Study'] == study
                    studies = plt.scatter(x[mask], y[mask], label=study_labels[i], color=colors[i], alpha=0.7, edgecolor='none', zorder=1)
                    studies_handles.append(study)

                #######NaOH
                x_NaOH = x1[0:len(naoh_columns)]

                #######Na2CO3
                # x_Na2CO3 = x1[len(naoh_columns):len(naoh_columns)*2]

                # x ranges for line plotting and confidence intervals
                xpred_l = linspace(x.min(), x.max(), 100)
                xm = x.mean()
                xs = x.std()
                min_val = x.min() - 0.5*xs
                max_val = max([x.max() + 0.5*xs, x1.max()+0.5*xs]) #min([xm + 3*xs, 1.25])
                if species in ['Duncanopsammia axifuga', 'Montastraea cavernosa']:
                    max_val = max([x.max() + 0.62*xs, x1.max()+0.62*xs]) 
                x1_smooth = np.linspace(min_val, max_val, 100)

                #Append # of datapoints, species names  and rate
                best_fits_studies.append(str(studies_handles).replace("'", "").strip("[]"))
                best_fits_n.append(len(x))
                best_fits_groups.append(rate_group['Group'].iloc[0])
                best_fits_species.append(species)
                                
                match = re.search(r'\[.*?\]', rate)
                try:
                    if match:
                        best_fits_rate.append(match.group(1))
                except IndexError:
                    best_fits_rate.append(rate)

                #Select model
                #Choose model based on p value
                smo = np.argmin([lp,xp])

                # Manually fix Halimeda opuntia fix since the fit is upside down threshold,
                # and linear is a better fit
                if species in ['Halimeda opuntia']:
                    smo=0

                if smo==0:
                    best_fits.append('linear')
                    best_fits_exponent.append(np.nan)
                    best_fits_slope.append(lparams[1])
                    best_fits_intercept.append(lparams[0])
                    best_fits_p.append(lp)
                    best_fits_r2.append(np.round(lrsqu, 2))
                    best_fits_rmse.append(np.round(lrmse, 4))

                    ypred_lin_l = fl(lparams, xpred_l)
                    plt.plot(xpred_l, ypred_lin_l, color='xkcd:grey')#, label='Experimental data species response')
                    ypred_smooth_l = fl(lparams, x1_smooth)
                    plt.plot(x1_smooth, ypred_smooth_l, color='xkcd:grey', linestyle='--')
                    y_NaOH = fl(lparams, x_NaOH)
                    # y_Na2CO3  = fl(lparams, x_Na2CO3)

                    # confidence interval
                    x1_smooth_const = sm.add_constant(x1_smooth)
                    conf_int = model_lin.get_prediction(x1_smooth_const).conf_int(obs=True, alpha=0.1)
                    upper_lin = conf_int[:, 1]
                    lower_lin = conf_int[:, 0]
                    plt.fill_between(x1_smooth, upper_lin, lower_lin, color='xkcd:grey', alpha=0.3, edgecolor='None', label = '90% Prediction Interval')

                    # pre-industrial rate
                    pre_industrial_rate = fl(lparams,df2_['TA-DIC pre-industrial'].values[0])
                    best_fits_pre_industrial.append(pre_industrial_rate)
                    best_fits_baseline_pi.append(df2_['TA-DIC pre-industrial'].values[0])
                    best_fits_TA_pi.append(df2_['TA pre-industrial'].values[0])
                    best_fits_DIC_pi.append(df2_['DIC pre-industrial'].values[0])

                    # confidence interval
                    # Combine x1_smooth, lower_lin, and upper_lin into a 2D array (each row is [x, lower, upper])
                    error_pi_data = np.column_stack((x1_smooth, lower_lin, upper_lin))
                    
                    # Find the index of the closest value in x1_smooth to the pre-industrial value
                    closest_index_pi = np.abs(x1_smooth - df2_['TA-DIC pre-industrial'].values[0]).argmin()
                    
                    # Extract the corresponding values
                    closest_x = x1_smooth[closest_index_pi]
                    closest_lower_pi = lower_lin[closest_index_pi]
                    closest_upper_pi = upper_lin[closest_index_pi]

                    # current rate
                    current_rate = fl(lparams,df2_['TA-DIC current'].values[0])
                    best_fits_current.append(current_rate)
                    best_fits_baseline.append(df2_['TA-DIC current'].values[0])
                    best_fits_TA_baseline.append(df2_['TA current'].values[0])
                    best_fits_DIC_baseline.append(df2_['DIC current'].values[0])

                    best_fits_NaOH_50_uneq.append(fl(lparams,df2_['For calc. TA-DIC +50 umol/kg NaOH uneq.'].values[0]))
                    best_fits_NaOH_50_eq.append(fl(lparams,df2_['For calc. TA-DIC +50 umol/kg NaOH eq.'].values[0]))

                    best_fits_Na2CO3_50_uneq.append(fl(lparams,df2_['For calc. TA-DIC +50 umol/kg Na2CO3 uneq.'].values[0]))
                    best_fits_Na2CO3_50_eq.append(fl(lparams,df2_['For calc. TA-DIC +50 umol/kg Na2CO3 eq.'].values[0]))
                                                
                    # Find the index of the closest value in x1_smooth to the pre-industrial value
                    closest_index = np.abs(x1_smooth - df2_['TA-DIC current'].values[0]).argmin()
                    
                    # Extract the corresponding values
                    closest_x = x1_smooth[closest_index]
                    closest_lower = lower_lin[closest_index]
                    closest_upper = upper_lin[closest_index]

                    best_fits_lower.append(closest_lower)
                    best_fits_upper.append(closest_upper)
                    best_fits_lower_pi.append(closest_lower_pi)
                    best_fits_upper_pi.append(closest_upper_pi)

                    # Compute ± range for current conditions
                    best_fits_range = (closest_upper - closest_lower) / 2
        
                    # Compute ± range for pre-industrial conditions
                    best_fits_range_pi = (closest_upper_pi - closest_lower_pi) / 2
                    
                    # Store results
                    best_fits_x_range.append(x.max()-x.min())

                elif smo==1:
                    best_fits.append('exponential')
                    best_fits_slope.append(xparam[1])
                    best_fits_intercept.append(xparam[0])
                    best_fits_exponent.append(xparam[2])
                    best_fits_p.append(xp)
                    best_fits_r2.append(np.round(xqsqu, 2))
                    best_fits_rmse.append(np.round(xrmse, 4))
                    
                    plt.plot(xpred_l, ypred_exp_l, color='xkcd:grey', label='Species experimental data regression')
                    ypred_smooth_l=fe(xparam, x1_smooth)
                    plt.plot(x1_smooth, ypred_smooth_l, color='xkcd:grey', linestyle='--')
                    y_NaOH = fe(xparam, x_NaOH)
                    # y_Na2CO3  = fe(xparam, x_Na2CO3)

                    ypred_l, ciLo, ciUp, p, rsqr = prediction_interval(f,xparam,x1_smooth,x.values,y.values,.1,True)
                    plt.fill_between(x1_smooth, ciLo, ciUp, color="xkcd:grey", alpha=0.2, edgecolor='None', label = '90% Prediction Interval')

                    # pre-industrial rate
                    pre_industrial_rate = fe(xparam, df2_['TA-DIC pre-industrial'].values[0])
                    best_fits_pre_industrial.append(pre_industrial_rate)
                    best_fits_baseline_pi.append(df2_['TA-DIC pre-industrial'].values[0])
                    best_fits_TA_pi.append(df2_['TA pre-industrial'].values[0])
                    best_fits_DIC_pi.append(df2_['DIC pre-industrial'].values[0])

                    # Combine x1_smooth, lower_lin, and upper_lin into a 2D array (each row is [x, lower, upper])
                    error_pi_data = np.column_stack((x1_smooth, ciLo, ciUp))
                    
                    # Find the index of the closest value in x1_smooth to the pre-industrial value
                    closest_index_pi = np.abs(x1_smooth - df2_['TA-DIC pre-industrial'].values[0]).argmin()
                    
                    # Extract the corresponding values
                    closest_x = x1_smooth[closest_index_pi]
                    closest_lower_pi = ciLo[closest_index_pi]
                    closest_upper_pi = ciUp[closest_index_pi]

                    # current rate
                    current_rate = fe(xparam, df2_['TA-DIC current'].values[0])
                    best_fits_current.append(current_rate)
                    best_fits_baseline.append(df2_['TA-DIC current'].values[0])
                    best_fits_TA_baseline.append(df2_['TA current'].values[0])
                    best_fits_DIC_baseline.append(df2_['DIC current'].values[0])

                    best_fits_NaOH_50_uneq.append(fe(xparam,df2_['For calc. TA-DIC +50 umol/kg NaOH uneq.'].values[0]))
                    best_fits_NaOH_50_eq.append(fe(xparam,df2_['For calc. TA-DIC +50 umol/kg NaOH eq.'].values[0]))

                    best_fits_Na2CO3_50_uneq.append(fe(xparam,df2_['For calc. TA-DIC +50 umol/kg Na2CO3 uneq.'].values[0]))
                    best_fits_Na2CO3_50_eq.append(fe(xparam,df2_['For calc. TA-DIC +50 umol/kg Na2CO3 eq.'].values[0]))
                                                
                    # Find the index of the closest value in x1_smooth to the pre-industrial value
                    closest_index = np.abs(x1_smooth - df2_['TA-DIC current'].values[0]).argmin()
                    
                    # Extract the corresponding values
                    closest_x = x1_smooth[closest_index]
                    closest_lower = ciLo[closest_index]
                    closest_upper = ciUp[closest_index]
 
                    best_fits_lower.append(closest_lower)
                    best_fits_upper.append(closest_upper)
                    best_fits_lower_pi.append(closest_lower_pi)
                    best_fits_upper_pi.append(closest_upper_pi)

                    # Compute ± range for current conditions
                    best_fits_range = (closest_upper - closest_lower) / 2
        
                    # Compute ± range for pre-industrial conditions
                    best_fits_range_pi = (closest_upper_pi - closest_lower_pi) / 2

                    # Store results
                    best_fits_x_range.append(x.max()-x.min())

                for i in range(len(x_NaOH)):  # Iterate through every other data point
                    plt.scatter(x_NaOH[i], y_NaOH[i],
                                color=cmap_NaOH[i],  # Use the sequential colormap value
                                label=f'{labels[i]} µmol/kg NaOH addition',
                                marker='*', edgecolor='grey', linewidth=0.5, s=50, zorder=10)
                    
                # for i in range(len(x_Na2CO3)):
                #     plt.scatter(x_Na2CO3[i], y_Na2CO3[i],
                #                 color=cmap_Na2CO3[i],
                #                 #label = f'{labels[i]} µmol/kg Na₂CO₃ addition',
                #                 marker='d', edgecolor='grey', linewidth=0.5, s=20, zorder=10, alpha=0.8)

                # find half calcification rate and then the Na2CO and NaOH rates
                yhalf = current_rate - current_rate*0.5
                yeps = current_rate*0.25

                a,b=x_NaOH.shape
                if b>0:
                    val1 = np.min(np.abs(y_NaOH-yhalf))
                    ind1 =  np.argmin(np.abs(y_NaOH-yhalf))
                    NaOH_threshold_species.append(species)
                    NaOH_threshold_groups.append(rate_group['Group'].iloc[0])
                    match = re.search(r'\[.*?\]', rate)
                    if match:
                        NaOH_threshold_rate.append(rate)
                    NaOH_thresholds.append(np.round(x_NaOH[ind1],3))

                # a,b=x_Na2CO3.shape
                # if b>0:
                #     val2 = np.min(np.abs(y_Na2CO3-yhalf))
                #     # if val2<=yeps:
                #     ind2 =  np.argmin(np.abs(y_Na2CO3-yhalf))
                #     Na2CO3_threshold_species.append(species)
                #     Na2CO3_threshold_groups.append(rate_group['Group'].iloc[0])
                #     match = re.search(r'\[.*?\]', rate)
                #     if match:
                #         Na2CO3_threshold_rate.append(rate)
                #     Na2CO3_thresholds.append(np.round(x_Na2CO3[ind2],3))
                #     # Na2CO3_threshold_amount.append(labels[ind2])
           
                # Function to find threshold values
                def find_thresholds(x_data, y_data, yhalf, yeps, threshold_species, thresholds, threshold_amount, species, labels):
                    if x_data.size > 0:
                        # Find indices where `y_data` is close to `yhalf`
                        indices = np.where(np.abs(y_data - yhalf) <= yeps)[0]
                        if len(indices) > 0:
                            for ind in indices:
                                NaOH_threshold_species.append(species)
                                NaOH_thresholds.append(np.round(x_data[ind], 3))
                                NaOH_threshold_amount.append(labels[ind])
                
                # Process NaOH data
                find_thresholds(x_NaOH, y_NaOH, yhalf, yeps, NaOH_threshold_species, NaOH_thresholds, NaOH_threshold_amount, species, labels)
                
                # Process Na2CO3 data
                # find_thresholds(x_Na2CO3, y_Na2CO3, yhalf, yeps, Na2CO3_threshold_species, Na2CO3_thresholds, Na2CO3_threshold_amount, species, labels)
 
                # Define intersection points
                x_current = df2_['TA-DIC current'].values[0]
                y_current = current_rate

                x_pre_industrial = df2_['TA-DIC pre-industrial'].values[0]
                y_pre_industrial = pre_industrial_rate
                
                # Get current axis limits to avoid zooming out
                if smo ==2:
                    y_std = np.std(y)  # Calculate the standard deviation of y
                    y_min = y.min() - (y_std * 0.5)
                    y_max = y.max() + (y_std * 0.5)
                    y_limits = (y_min, y_max)
                else:
                    y_limits = plt.gca().get_ylim()
                

                # Line width
                lw = 1.5  # Set this to your preferred line width
                
                # Plot horizontal line for current rate, stopping at x_current
                plt.plot([min_val, x_current], [y_current, y_current], color='xkcd:light orange', zorder=1.5, linestyle='-', linewidth=1.7, label='Current calcification rate')
                
                # Plot vertical line for Baseline current, stopping at y_current
                plt.plot([x_current, x_current], [y_limits[0], y_current], color='xkcd:light orange', zorder=1.5, linestyle='-', linewidth=1.7)

                # Plot horizontal line for pre-industrial rate, stopping at x_pre_industrial
                plt.plot([min_val, x_pre_industrial], [y_pre_industrial, y_pre_industrial], color='xkcd:deep pink',  zorder=1, linestyle='-', linewidth=1.5, label = 'Pre-industrial calcification rate')
                
                # Plot vertical line for Baseline pre-industrial, stopping at y_pre_industrial
                plt.plot([x_pre_industrial, x_pre_industrial], [y_limits[0], y_pre_industrial], color='xkcd:deep pink', zorder=1, linestyle='-', linewidth=1.5)
                
                plt.ylim(y_limits)      
                
                
            else:
                smo = 3
                studies_handles = []
                for i, study in enumerate(rate_group['Study'][~nan_mask].unique()):
                    mask = rate_group['Study'] == study
                    studies = plt.scatter(x[mask], y[mask], label=study, color=colors[i], alpha=0.7, edgecolor='none', zorder=0)
                    studies_handles.append(study)
                best_fits_studies.append(str(studies_handles).replace("'", "").strip("[]"))
                best_fits_n.append(len(x))
                best_fits_groups.append(rate_group['Group'].iloc[0])
                best_fits_species.append(species)
                match = re.search(r'\[.*?\]', rate)
                if match:
                    best_fits_rate.append(rate)

                best_fits.append('neutral')
                best_fits_p.append(np.nan)
                best_fits_slope.append('-')
                best_fits_intercept.append(np.nan)
                best_fits_exponent.append(np.nan)
                best_fits_r2.append(np.nan)
                best_fits_rmse.append(np.nan)
                best_fits_pre_industrial.append(np.nan)
                best_fits_current.append(np.nan)
                                
                best_fits_baseline.append(np.nan)
                best_fits_baseline_pi.append(np.nan)
                
                best_fits_lower.append(np.nan)
                best_fits_upper.append(np.nan)
                best_fits_lower_pi.append(np.nan)
                best_fits_upper_pi.append(np.nan)
                best_fits_x_range.append(np.nan)
                                                                               
                best_fits_TA_baseline.append(np.nan)
                best_fits_DIC_baseline.append(np.nan)
                best_fits_TA_pi.append(df2_['TA pre-industrial'].values[0])
                best_fits_DIC_pi.append(df2_['DIC pre-industrial'].values[0])

                best_fits_NaOH_50_uneq.append(np.nan)
                best_fits_NaOH_50_eq.append(np.nan)

                best_fits_Na2CO3_50_uneq.append(np.nan)
                best_fits_Na2CO3_50_eq.append(np.nan)
                                        
            plt.xlim([min_val, max_val])

            #Formatting
            species_formatted = species.replace(" ", "\u00A0")

            #Axes labels        
            plt.xlabel('TA-DIC', fontsize=9)
            
            if rate == rate1:
                plt.ylabel('Calc. rate ${CaCO}_3$ $[{mmol}$ ${g}^{-1}$ ${hr}^{-1}$]', fontsize=9)
                rate_label = 'Rate 1'
            if rate == rate2:
                plt.ylabel('Calc. rate ${CaCO}_3$ $[{mmol}$ ${m}^{-2}$ ${hr}^{-1}$]', fontsize=9)
                rate_label = 'Rate 2'
            if rate == rate3:
                plt.ylabel('Calc. rate ${CaCO}_3$ $[{mmol}$ ${m}^{-3}$ ${hr}^{-1}$]', fontsize=9)
                rate_label = 'Rate 3'
            if rate == rate4:
                plt.ylabel('Calc. rate $CaCO_3$ $[{mmol}$ $\#^{-1}$ $hr^{-1}]$', fontsize=9)
                rate_label = 'Rate 4'
            if rate == rate5:
                plt.ylabel('Calc. rate ${CaCO}_3$ $[{mmol}$ ${hr}^{-1}$]', fontsize=9)
                rate_label = 'Rate 5'
            if rate == rate6:
                plt.ylabel('Calc. rate ${CaCO}_3$ $[{mmol}$ ${cm}^{-2}$]', fontsize=9)
                rate_label = 'Rate 6'
            if rate == rate7:
                plt.ylabel('Calc. rate ${CaCO}_3$ $[{\%}$ ${hr}^{-1}$]', fontsize=9)
                rate_label = 'Rate 7'

            # Define a formatter function for the y-axis
            def scientific_format(y, pos):
                return f"{y:.1E}"  # Customize the number of decimal points here (e.g., .2 for two decimal points)
            
            # Apply the formatter with plt.gca() to ensure uniform scientific formatting
            plt.gca().yaxis.set_major_formatter(ticker.FuncFormatter(scientific_format))

            #Formatting
            species_formatted = species.replace(" ", "\u00A0")
            plt.title('{}\n${}$'.format(rate_group['Group'].iloc[0], species_formatted),fontsize=9)

            # #Axes labels        
            plt.xlabel('TA-DIC', fontsize=9)
            
            #Complete legend
            plt.tight_layout()
            
            if (smo == 1) or (smo == 0 and lparams[1] > 0) and (smo != 3):
                plt.savefig(filepath+'/Figures/{} ({}) {}'.format(rate_group['Group'].iloc[0], species_formatted, rate_label)+'.png', bbox_inches='tight')
                plt.show()

#%%
#Make df with responses and stats
best_fits_df = pd.DataFrame({
    'Studies': best_fits_studies,
    '# datapoints': best_fits_n,
    'Group': best_fits_groups,
    'Species': best_fits_species,
    'Rate unit': best_fits_rate,
    'Response': best_fits,
    'Slope': best_fits_slope,
    'Intercept': best_fits_intercept,
    'Exponent': best_fits_exponent,
    'p-value': best_fits_p,
    'R2': best_fits_r2,
    'RMSE': best_fits_rmse,
    'Current TA': best_fits_TA_baseline,
    'Current DIC': best_fits_DIC_baseline,
    'TA-DIC current': best_fits_baseline,
    'Range TA-DIC current': best_fits_x_range,
    'Current rate': best_fits_current,
    'Pre-industrial TA': best_fits_TA_pi,
    'Pre-industrial DIC': best_fits_DIC_pi,
    'TA-DIC pre-industrial': best_fits_baseline_pi,
    'Pre-industrial rate': best_fits_pre_industrial,
    'Calc. rate +50 NaOH uneq': best_fits_NaOH_50_uneq,
    'Calc. rate +50 NaOH eq': best_fits_NaOH_50_eq,
    'Calc. rate +50 Na2CO3 uneq': best_fits_Na2CO3_50_uneq,
    'Calc. rate +50 Na2CO3 eq': best_fits_Na2CO3_50_eq
    })

best_fits_df['Slope'] = pd.to_numeric(best_fits_df['Slope'], errors='coerce')
best_fits_df = best_fits_df[(best_fits_df['Response'] != 'neutral')]
best_fits_df = best_fits_df[~((best_fits_df['Response'] == 'linear') & (best_fits_df['Slope'] < 0))]
best_fits_df = best_fits_df.sort_values(by=['Group', 'Species'], ascending=[True, True])

