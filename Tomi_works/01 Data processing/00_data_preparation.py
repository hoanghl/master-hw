

#### Mini-project #####

# Make sure that .csv files are in your working directory!

## Load needed packages. ##
import pandas as pd

## Load datasets first each to its own data frame. ##
df_rp = pd.read_csv('global-data-on-sustainable-energy.csv')
print(df_rp.head(10))
df_tp = pd.read_csv('eleccap.csv')
#print(df_tp.head(10))
df_mp = pd.read_csv('pubfin.csv')
#print(df_mp.head(10))


## Merge data frames into a single data frame. ##

df_merged = df_tp.merge(df_mp, how = 'outer')
#print(df_merged.head(10))
df_merged = df_rp.merge(df_merged, how = 'outer')
#print(df_merged.head(10))


# Let's take only those columns which we need for the final data frame.

columns = ['Entity', 'Year', 'gdp_growth', 'gdp_per_capita', 'Electricity from renewables (TWh)', 'Capacity(MW)', 'Investment(USDmn)']
df = df_merged[columns]
#print(df.head(10))

# Let's export the original data frame into a .csv-file.
# The purpose is to check how our data imputation worked. 

df_org = df.copy()

# Then we can limit our data to years 2011-2020.
df_org = df_org[df_org['Year'] >= 2011]
df_org = df_org[df_org['Year'] <= 2020]

df_org.to_csv('data_merged.csv', index = False)

# Let's check how many NaN-values we have in the 3 interesting columns.

rp_nan_count = df['Electricity from renewables (TWh)'].isna().sum()
print(f'number of missing values of RP before {rp_nan_count}')

tp_nan_count = df['Capacity(MW)'].isna().sum()
print(f'number of missing values of TP before {tp_nan_count}')

mp_nan_count = df['Investment(USDmn)'].isna().sum()
print(f'number of missing values of MP before {mp_nan_count}')

print("")

# Let's use interpolation and forward fill to replace missing values.

print(f'Interpolation')

df_int = df.copy()
df_int.interpolate() # Interpolate first
df_int.fillna(method = 'ffill', inplace = True) # If still missing values, fill them with forwar fill

# Then we can limit our data to years 2011-2020.
df_int = df_int[df_int['Year'] >= 2011]
df_int = df_int[df_int['Year'] <= 2020]

rp_nan_count = df_int['Electricity from renewables (TWh)'].isna().sum()
print(f'number of missing values of RP after {rp_nan_count}')

tp_nan_count = df_int['Capacity(MW)'].isna().sum()
print(f'number of missing values of TP after {tp_nan_count}')

mp_nan_count = df_int['Investment(USDmn)'].isna().sum()
print(f'number of missing values of MP after {mp_nan_count}')

# Then we can export the data frame into a .csv-file 

df_int.to_csv('data_merged_int.csv', index = False)

"""
# Let's use forward fill to replace missing values.

print(f'Forward fill')

df_ffill = df.copy()
df_ffill.fillna(method = 'ffill', inplace = True)

# Then we can limit our data to years 2011-2020.
df_ffill = df_ffill[df_ffill['Year'] >= 2011]
df_ffill = df_ffill[df_ffill['Year'] <= 2020]

rp_nan_count = df_ffill['Electricity from renewables (TWh)'].isna().sum()
print(f'number of missing values of RP after {rp_nan_count}')

tp_nan_count = df_ffill['Capacity(MW)'].isna().sum()
print(f'number of missing values of TP after {tp_nan_count}')

mp_nan_count = df_ffill['Investment(USDmn)'].isna().sum()
print(f'number of missing values of MP after {mp_nan_count}')

# Then we can export the data frame into a .csv-file 

df_ffill.to_csv('data_merged_ffill.csv', index = False)

print('')

"""

