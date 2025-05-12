##########################################################
## Compute current and pre-industrial conditions, conceptually
## add alkalinity (in the form of NaOH and Na2CO3) 
##########################################################
## Author: H. van de Mortel
## Version: 1.0
## Maintainer: H. van de Mortel
## Email: vandemortel.hanna@gmail.com
##########################################################

import pandas as pd
import PyCO2SYS as pyco2

# Can consult PyCO2SYS manual:
# https://pyco2sys.readthedocs.io/en/latest/co2sys_nd/


filepath = "C:/Users/hanna/Documents/GitHub/OAE_mitigation_of_OA/"

# Import data:
# df1 should include (at minimum) two columns with carb parameters to solve 
# the carbonate system (in this case TA and DIC), a column with species names 
# and a column with calcification rate data. If species and/or rate units
# are all uniform, the for-loop can be simplified.

df1 = pd.read_excel(filepath + "master_data.xlsx", 
                    sheet_name='Sheet1')

# df_c includes the same columns as df1, but with the experimental control 
# conditions of each study (written out per species).

df_c = pd.read_excel(filepath + "control_conditions.xlsx")

#%% Variables and grouping

# Define different rate units, so they can be handled separately in for-loop
# (These are all separate column headers in df1; change accordingly)
rate1 = 'Calc rate [mmol/g/hr] (compiled)'
rate2 = 'Calc rate [mmol/m**2/hr] (compiled)'
rate3 = 'Calc rate [mmol/m**3/hr] (compiled)'
rate4 = 'Calc rate [mmol/#/hr] (compiled)'
rate5 = 'Calc rate [mmol/hr] (compiled)'
rate6 = 'Calc rate [mmol/cm**2] (compiled)'
rate7 = 'Calc rate [%/hr] (compiled)'

# Create list for for-loop
rates = [rate1, rate2, rate3, rate4, rate5, rate6, rate7] 

# Group by species for for-loop
grouped = df1.groupby('Species (compiled)')

# Estimation silicate and phosphate concentrations (minor influence)
silicate = 50  # umol/kg (= 30 uM)
phosphate = 0.5  # umol/kg (= 0.05 mg/L)

#%% TA ADDITION

# Create an empty list to store DataFrames
results_list = []

# Loop through species
for species, group in grouped:
    # Loop through different rate units
    for rate in rates:          
        if rate in group.columns and not group[rate].dropna().empty:
            
            # Remove nan values and make sure data is all numeric
            rate_group = group[group[rate].notna()]
            rate_group.loc[:, rate] = pd.to_numeric(rate_group[rate], errors='coerce')

            # Filter relative experimental control conditions for this species
            species_ctrl = df_c[df_c['Species'] == species]

### CURRENT & PRE-INDUSTRIAL CONDITIONS ###
            # Define parameters for PyCO2SYS, using experimental control sal and temp
            parms1 = dict(
                         salinity = species_ctrl['Sal ctrl'].values[0], 
                         temperature = species_ctrl['Temp ctrl'].values[0],
                         total_silicate=silicate, total_phosphate=phosphate,
                         opt_pH_scale=1, #total scale
                         opt_k_carbonic=10, #Lueker 2000
                         opt_total_borate=2 #Lee 2010
                         #Sulfuric acid dissociation = Dickson 1990
                         )
            

            #Compute current conditions baselines using experimental control TA and DIC
            baseline = pyco2.sys(**parms1, 
                                par1 = species_ctrl['TA ctrl'].values[0], par1_type = 1,
                                par2 = species_ctrl['DIC ctrl'].values[0], par2_type = 2,
)
            
            # Store current conditions values as variables
            TA_baseline = baseline['alkalinity']
            DIC_baseline = baseline['dic']
            pH_baseline = baseline['pH']
            pCO2_baseline = baseline['pCO2']
            
            # Compute pre-industrial conditions using 
            # current TA and pCO2 of current - 142 ppm
            pi_conditions = pyco2.sys(**parms1, 
                                par1 = TA_baseline, par1_type=1,
                                par2 = pCO2_baseline-142, par2_type = 4,
                                )
            
            # Store pre-industrial values as variables
            TA_pi = pi_conditions['alkalinity']
            DIC_pi = pi_conditions['dic']
            pH_pi = pi_conditions['pH']
            pCO2_pi = pi_conditions['pCO2']
            
            # Save variables in results_dict
            results_dict1 = {
            'Group': rate_group['Group'].iloc[0],
            'Species': species,
            'Rate': rate,
            'TA-DIC current': TA_baseline-DIC_baseline,
            'TA current': TA_baseline,
            'DIC current': DIC_baseline,
            'pH current': pH_baseline,
            'pCO2 current': pCO2_baseline,
            'TA-DIC pre-industrial': TA_pi-DIC_pi,
            'TA pre-industrial': TA_pi,
            'DIC pre-industrial': DIC_pi,
            'pH pre-industrial': pH_pi,
            'pCO2 pre-industrial': pCO2_pi
            }
            
###UNEQUILIBRATED###
#(for plotting)
#NaOH
                
            for i in list(range(50, 201, 50)):
                parms2 = dict(#salinity=rate_group['Sal'].mean(), temperature = rate_group['Temp [°C]'].mean(),
                        salinity = species_ctrl['Sal ctrl'].values[0], temperature = species_ctrl['Temp ctrl'].values[0],
                        total_silicate=silicate, total_phosphate=phosphate,
                          opt_pH_scale=1, #total scale
                          opt_k_carbonic=10, #Lueker 2000
                          opt_total_borate=2 #Lee 2010
                          #Sulfuric acid dissociation = Dickson 1990
                          )
                
                #########NaOH#########
                # Compute carbonate system for each addition of NaOH
                results_NaOH = pyco2.sys(**parms2, 
                                    par1 = TA_baseline+i, par1_type = 1,
                                    par2 = DIC_baseline, par2_type = 2)

                # Populate results_dict with carb system variables of interest
                # Define carb system variables of interest and corresponding 
                # result fields
                fields1 = [
                    ('TA-DIC', 'alkalinity', lambda r: r['alkalinity'] - r['dic']),
                ]

                # Populate results_dict using a loop
                for field_suffix, result_key, result_extractor in fields1:
                    column_name = f'NaOH {field_suffix} +{i:.0f} umol/kg'
                    results_dict1[column_name] = result_extractor(results_NaOH)

###UNEQUILIBRATED###
#(for exact TA addition to reach pre-industrial)
#NaOH & Na2CO3

            dNaOH_uneq = (TA_pi - DIC_pi) - (TA_baseline - DIC_baseline)
            dNa2CO3_uneq = 2 * dNaOH_uneq
            results_dict1.update({
                'dNaOH_uneq': dNaOH_uneq,
                'dNa2CO3_uneq': dNa2CO3_uneq
            })

###EQUILIBRATED###
#(for exact TA addition to reach pre-industrial)

            Alk_pi = TA_pi-DIC_pi               # Pre-industrial Alk*
            Alk_t0 = TA_baseline-DIC_baseline   # Initial Alk*
            CDReff = 0.8                        # CDR efficiency
            etamax = 0.832                       # Initial ηmax
 
            if species in ['Limacina helicina']:
                etamax = 0.904
                
            dNaOH_eq = (Alk_pi-Alk_t0)/(1-CDReff*etamax)
            dNa2CO3_eq = (Alk_pi-Alk_t0)/(0.5-CDReff*(etamax-0.5))

###Calc. rate for +50 NaOH uneq. and eq.###
            #For NaOH
            z = 2.99
            if species == 'Limacina helicina':
                z= 3.61
            
            i_values = [50, 50/z]
            i_labels = ['+50 umol/kg NaOH uneq.', '+50 umol/kg NaOH eq.']

            for i, label in zip (i_values, i_labels):
                parms2 = dict(
                         salinity = species_ctrl['Sal ctrl'].values[0], temperature = species_ctrl['Temp ctrl'].values[0],
                         total_silicate=silicate, total_phosphate=phosphate,
                         opt_pH_scale=1, #total scale
                         opt_k_carbonic=10, #Lueker 2000
                         opt_total_borate=2 #Lee 2010
                         #Sulfuric acid dissociation = Dickson 1990
                         )
                
                #########NaOH#########
                # Compute carbonate system for each addition of NaOH
                results_NaOH = pyco2.sys(**parms2, 
                                    par1 = TA_baseline+i, par1_type = 1,
                                    par2 = DIC_baseline, par2_type = 2)

                # Define carb system variables of interest and corresponding 
                # result fields
                fields1 = [
                    ('TA-DIC', 'alkalinity', lambda r: r['alkalinity'] - r['dic'])]

                # Populate results_dict using a loop
                for field_suffix, result_key, result_extractor in fields1:
                    column_name = f'For calc. {field_suffix} {label}'
                    results_dict1[column_name] = result_extractor(results_NaOH)
                    
###Calc. rate for +0 Na2CO3 uneq. and eq.###
            #For Na2CO3
            z = 2.13

            if species == 'Limacina helicina':
                z= 2.83
            
            i_values = [50, 50/z]
            i_labels = ['+50 umol/kg Na2CO3 uneq.', '+50 umol/kg Na2CO3 eq.']

            for i, label in zip (i_values, i_labels):
                parms2 = dict(
                         salinity = species_ctrl['Sal ctrl'].values[0], temperature = species_ctrl['Temp ctrl'].values[0],
                         total_silicate=silicate, total_phosphate=phosphate,
                         opt_pH_scale=1, #total scale
                         opt_k_carbonic=10, #Lueker 2000
                         opt_total_borate=2 #Lee 2010
                         #Sulfuric acid dissociation = Dickson 1990
                         )
                #########NaOH#########
                # Compute carbonate system for each addition of NaOH
                results_Na2CO3 = pyco2.sys(**parms2, 
                                    par1 = TA_baseline+i, par1_type = 1,
                                    par2 = DIC_baseline+(0.5*i), par2_type = 2)

                # Define carb system variables of interest and corresponding 
                # result fields
                fields1 = [
                    ('TA-DIC', 'alkalinity', lambda r: r['alkalinity'] - r['dic'])]

                # Populate results_dict using a loop
                for field_suffix, result_key, result_extractor in fields1:
                    column_name = f'For calc. {field_suffix} {label}'
                    results_dict1[column_name] = result_extractor(results_Na2CO3)
                    
            results_dict1.update({
                'dNaOH_eq': dNaOH_eq,
                'dNa2CO3_eq': dNa2CO3_eq,
                'etamax': etamax,
                'etamax and CDReff TA-DIC line': (TA_baseline+dNaOH_uneq)- (DIC_baseline+(0.8*0.85*dNaOH_uneq)),
                'Conversion factor dNaOH': dNaOH_eq/dNaOH_uneq,
                'Conversion factor dNa2CO3': dNa2CO3_eq/dNa2CO3_uneq,
                'Max TA': rate_group['AT [µmol/kg] (compiled)'].max(),
                'Min TA': rate_group['AT [µmol/kg] (compiled)'].min(),
                'Median TA': rate_group['AT [µmol/kg] (compiled)'].median()

            })            

            # Append results into results_dict
            results_list.append(results_dict1)  
            
# Concatenate all DataFrames in the list into a single DataFrame
final_results = pd.DataFrame(results_list)
final_results = final_results.sort_values(by=['Group', 'Species'], ascending=[True, True])

# Save df to Excel file
output_file_path = filepath+"/output_TA_addition.xlsx"
final_results.to_excel(output_file_path, index=False)