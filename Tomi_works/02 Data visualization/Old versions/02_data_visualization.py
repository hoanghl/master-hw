

##### Mini-project #####

# Change working directory such that data files are in the workind directory.

## Load needed packages. ##
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pycountry_convert as pc

def country_to_continent(country_name):
    try:
        country_alpha2 = pc.country_name_to_country_alpha2(country_name)
        country_continent_code = pc.country_alpha2_to_continent_code(country_alpha2)
        country_continent_name = pc.convert_continent_code_to_continent_name(country_continent_code)
        return country_continent_name
    except:
        return None

## Load processed dataset. ##
df = pd.read_csv('data_merged_int.csv')
#print(df.head(10))
print

# Let's determine first the continent of each country.

n = len(df)
continents = []

for i in range(n):
    country = df['Entity'].iloc[i]
    #print(country)
    continent = country_to_continent(country)
    continents.append(continent)
#print(continents)

df['Continent'] = continents
#print(df.head(10))

# Lets first calculate the change in three variables we are interested in.

country = df['Entity'].iloc[0]
rp = [0]
tp = [0]
mp = [0]

n = len(df)
#print(n)
#print(df.iloc[5])

for i in range(n - 1):
    if df['Entity'].iloc[i + 1] != country:
        rp.append(np.nan)
        tp.append(np.nan)
        mp.append(np.nan)
    else:
        rp_im1 = df['Electricity from renewables (TWh)'].iloc[i]
        rp_i = df['Electricity from renewables (TWh)'].iloc[i + 1]
        tp_im1 = df['Capacity(MW)'].iloc[i]
        tp_i = df['Capacity(MW)'].iloc[i + 1]
        mp_im1 = df['Investment(USDmn)'].iloc[i]
        mp_i = df['Investment(USDmn)'].iloc[i + 1]
        if rp_im1 == 0:
            rp.append(np.inf)
        else:
            rp.append(rp_i/rp_im1)
        if tp_im1 == 0:
            tp.append(np.inf)
        else:
            tp.append(tp_i/tp_im1)
        if mp_im1 == 0:
            mp.append(np.inf)
        else:
            mp.append(mp_i/mp_im1)

    country = df['Entity'].iloc[i + 1]
    

df['RP'] = rp
df['TP'] = tp
df['MP'] = mp

# If there is a case of "divede-by-zero", lets replace these values with 0.
df = df.replace([np.inf], 0)

rei = []

# Let's calculate the value of the REI-index for all countries in the sample.

for i in range(n):
    val = (df['RP'].iloc[i] + df['TP'].iloc[i] + df['MP'].iloc[i])/3
    rei.append(val)

df['REI'] = rei

#print(df.head(20))
#print(df[df['Entity'] == 'United States of America'])

sample = ['Australia', 'Canada', 'Finland', 'France', 'Germany']
#Ireland, Israel, Italy, Japan, Norway, Portugal, South Africa, Spain, Sweden, United Kingdom, United States, Republic of Korea

df_sample = df.copy()
df_sample = df_sample[df_sample['Entity'].isin(sample)]
#print(df_sample.head(10))

# Let's plot now graphs.

for country in sample:
    df_fig  = df_sample[df_sample['Entity'] == country]
    #print(df_fig)
    df_fig.plot('Year', y = ['RP', 'TP', 'MP', 'REI'], kind = 'bar')
    plt.title(country)
    name = 'timeseries_' + country + '.png'
    plt.savefig(name)

df_pivot = df_sample.pivot(index = "Year", columns = "Entity", values = "REI")
print(df_pivot.head(10))

#df_pivot.plot(kind = 'bar')
#plt.title('REI for different countries')

#plt.show()

##### Take only the latest year for the analysis. #####

df_latest = df.copy()
df_latest = df_latest[df_latest['Year'] == 2020]
print(df_latest.head(10))

df_latest = df_latest.sort_values(by = ['Continent', 'REI'], ascending = [True, False])
#print(df_latest[df['Continent'] == 'Europe'].head(10))

# Plot the REI of the 10 highest from each continent.

continents_unique = list(df.Continent.unique())
print(continents_unique)

#print(df_latest[df['Continent'] == 'Europe'].iloc[:10])
#df_latest[df['Continent'] == 'Europe'].iloc[:10].plot('Entity', y = ['REI'], kind = 'bar')
#plt.title('Europe')
#plt.show()

for cont in continents_unique:
    if cont != None:
        df_latest[df['Continent'] == cont].iloc[:10].plot('Entity', y = ['RP', 'TP', 'MP', 'REI'], kind = 'bar', figsize=(10,6))
        plt.title(f'10 highest REI in {cont} in 2020')
        name = 'higest_REI_' + cont + '.png'
        plt.savefig(name)

#plt.show()
