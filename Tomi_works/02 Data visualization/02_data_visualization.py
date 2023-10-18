

##### Mini-project #####

# Change working directory such that data files are in the workind directory.

## Load the needed packages. ##
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pycountry_convert as pc

#### 1. Pre-prosessing of the input data. ####

#######################################################################

# Function which finds the continent for the country.
def country_to_continent(country_name):
    try:
        country_alpha2 = pc.country_name_to_country_alpha2(country_name)
        country_continent_code = pc.country_alpha2_to_continent_code(country_alpha2)
        country_continent_name = pc.convert_continent_code_to_continent_name(country_continent_code)
        return country_continent_name
    except:
        return None

# Function which returns 2-letter country code for the country.
def country_alpha2(countries):
    try:
        alpha2_codes = []
        for country in countries:
            alpha2 = pc.country_name_to_country_alpha2(country)
            alpha2_codes.append(alpha2)
        return alpha2_codes
    except:
        return None


#######################################################################

# Load the processed dataset.
df = pd.read_csv('data_merged_int.csv')
#print(df.head(10))

#######################################################################

# Let's find out the continent of each country.
n = len(df)
continents = []

for i in range(n):
    country = df['Entity'].iloc[i]
    #print(country)
    continent = country_to_continent(country)
    continents.append(continent)
#print(continents)

df['Continent'] = continents
#print(df.Continent.unique())

#### UNCOMMENT IF THERE IS A NEED TO CHECK WHICH COUNTRIES HAVE NO CONTINENT. ####
# Let's check which entities do not belong to any continent.
#print('These entities do not beling to any continent:')
#print(df[df['Continent'].isnull()].Entity.unique())
#print('\n')
#print(df.head(10))

#######################################################################

# Lets now calculate the change in three variables we are interested in.
country = df['Entity'].iloc[0]
rp = [0]
tp = [0]
mp = [0]

n = len(df)
#print(n)
#print(df.iloc[5])

for i in range(n - 1):
    # Since we calculate the differences, the first date does not have any value since there is nothing to compare against.
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
        # This is just a check if we end up dividing with zero then the ratio is infinite.
        # We will deal with this later on.
        if rp_im1 == 0:
            rp.append(np.inf)
        # If not, then we can calculate the ratio.
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
    
# Then we add the calculated ratios to the dataframe.

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

#######################################################################

#################### 2. Visualization of the data. ####################

########################### 2.1 Overall. ##############################

## Plot the evolution of all three variables. ##

# Let's first create a copy of the dataframe.
df_agg = df.groupby('Year').sum(numeric_only = True)
#print(df_agg.head())

#### UNCOMMENT WHEN YOU WANT TO SAVE OR VIEW PLOTS. ####
df_agg.plot(y = ['Electricity from renewables (TWh)'], kind = 'bar')
plt.savefig('overall_elec_timeseries.png')
df_agg.plot(y = ['Capacity(MW)'], kind = 'bar')
plt.savefig('overall_cap_timeseries.png')
df_agg.plot(y = ['Investment(USDmn)'], kind = 'bar')
plt.savefig('overall_inv_timeseries.png')
#plt.show()

#######################################################################

## Plot the distribution of the three variables in different continents in 2020. ##

# Lets crate a copy of the dataframe.
df_all = df[df['Year'] == 2020]

df_all_agg = df_all.groupby('Continent').sum(numeric_only = True)
#print(df_all.head(11))

#### UNCOMMENT WHEN YOU WANT TO SAVE OR VIEW PLOTS. ####
df_all_agg.plot.pie(y = 'Electricity from renewables (TWh)', startangle = 90, counterclock = False, title = 'Year 2020')
plt.savefig('overall_elec_pie_2020.png')
df_all_agg.plot.pie(y = 'Capacity(MW)', startangle = 90, counterclock = False, title = 'Year 2020')
plt.savefig('overall_cap_pie_2020.png')
df_all_agg.plot.pie(y = 'Investment(USDmn)', startangle = 90, counterclock = False, title = 'Year 2020')
plt.savefig('overall_inv_pie_2020.png')
#plt.show()

#######################################################################

##### 10 highest REI from each continent in 2020. #####

#df_all = df_all.reset_index()
#print(df_all.head(10))

df_all = df_all.sort_values(by = ['Continent', 'REI'], ascending = [True, False])
#print(df_all[df['Continent'] == 'Europe'].head(10))

# Plot the REI of the 10 highest from each continent.

continents_unique = list(df.Continent.unique())
#print(continents_unique)

#### UNCOMMENT WHEN YOU WANT TO SEE THE GRAPH! ####
for cont in continents_unique:
    if cont != None:
        df_fig = df_all.copy()
        df_fig = df_fig[df_fig['Continent'] == cont].iloc[:10]
        #df_fig = df_fig.reset_index()
        #print(df_fig.head(10))
        df_fig.plot('Entity', y = ['RP', 'TP', 'MP', 'REI'], kind = 'bar')
        plt.title(f'10 highest REI in {cont} in 2020')
        name = cont + '_highest_REI' + '.png'
        plt.savefig(name)

#plt.show()

###################### 2.2 Continent-by-continent #####################

continents = continents_unique[: len(continents_unique) - 1]
#print(continents)

for i in range(len(continents)):
    continent = continents[i]

    ############################### Continent ################################

    ## Find out which countries in each continent have the higest REI in 2020. ##

    # Lets crate a copy of the dataframe.
    df_eu = df[df['Continent'] == continent]
    df_eu_2020 = df_eu[df_eu['Year'] == 2020]

    # Lets sort values based on the REI.
    df_eu_2020 = df_eu_2020.sort_values(by = ['REI'], ascending = [False])

    # 10 counties in Europe with highest REI in 2020.
    eu_top_10 = list(df_eu_2020['Entity'].iloc[:10]) 
    #print(eu_top_10)

    # 2-letter codes for those countries.
    eu_top_10_a2 = country_alpha2(eu_top_10)
    #print(eu_top_10_a2)

    #######################################################################

    # Plot the evolution of the three variables for the countries with 
    # the highest (top-10) REI in each continent in 2020.

    df_eu_top10 = df_eu[df_eu['Entity'].isin(eu_top_10)]
    df_eu_top10 = df_eu_top10.reset_index()
    #df_eu_top10 = df_eu_top10[df_eu_top10['Entity'] in eu_top_10]
    #print(df_eu_top10.head(10))

    # Let's plot now graphs.
    for country in eu_top_10:
        df_fig  = df_eu_top10[df_eu_top10['Entity'] == country]
        df_fig = df_fig.reset_index()
        #print(df_fig)
        #### UNCOMMENT WHEN YOU WANT TO SEE THE GRAPH! ####
        df_fig.plot('Year', y = ['RP', 'TP', 'MP', 'REI'], kind = 'bar')
        plt.title(country)
        name = continent + "_" + country + '_timeseries' + '.png'
        plt.savefig(name)
        #plt.show()
        
    df_pivot = df_eu_top10.pivot(index = "Year", columns = "Entity", values = "REI")
    #print(df_pivot.head(10))

    #### UNCOMMENT WHEN YOU WANT TO SEE THE GRAPH! ####
    df_pivot.plot(kind = 'bar')
    plt.title(f'REI for different countries in {continent}')
    name = continent + "_" + "top10" + '_timeseries' + '.png'
    plt.savefig(name)
    #plt.show()
    
    #######################################################################

    ##### This block is not needed!!!! #####

    eu_names_all = []

    # Lets rename countries which are not in top-10 for the whole dataframe.
    for i in range(len(df_eu)):
        name = df_eu['Entity'].iloc[i]
        if name not in eu_top_10:
            eu_names_all.append('Other')
        else:
            eu_names_all.append(name)

    df_eu = df_eu.reset_index()
    df_eu['Entity t10'] = eu_names_all

    #######################################################################

    # Plot the distribution of the three variables in Europe in 2020
    # in a pie chart where top-10 countries have their name and the rest
    # are groudes as 'Other'.

    eu_names_2020 = []

    # Lets rename countries which are not in top-10 for the 2020 dataframe.
    for i in range(len(df_eu_2020)):
        name = df_eu_2020['Entity'].iloc[i]
        if name not in eu_top_10:
            eu_names_2020.append('Other')
        else:
            alpha2_name = pc.country_name_to_country_alpha2(name)
            eu_names_2020.append(alpha2_name)

    df_eu_2020 = df_eu_2020.reset_index()
    df_eu_2020['Entity t10'] = eu_names_2020

    df_eu_2020_agg = df_eu_2020.groupby('Entity t10').sum(numeric_only = True)
    #print(df_eu_2020_agg.head(11))

    #### UNCOMMENT WHEN YOU WANT TO SEE/SAVE THE GRAPH! ####
    df_eu_2020_agg.plot.pie(y = 'Electricity from renewables (TWh)', startangle = 90, counterclock = False, title = continent + ' Year 2020')
    plt.title(f'Electricity from renewables (TWh) in {continent}')
    name = continent + "_elec_pie_2020.png"
    plt.savefig(name)
    
    df_eu_2020_agg.plot.pie(y = 'Capacity(MW)', startangle = 90, counterclock = False, title = 'Year 2020')
    plt.title(f'Capacity (MW) in {continent}')
    name = continent + "_cap_pie_2020.png"
    plt.savefig(name)
    
    df_eu_2020_agg.plot.pie(y = 'Investment(USDmn)', startangle = 90, counterclock = False, title = 'Year 2020')
    plt.title(f'Investments (USDmn) in {continent}')
    name = continent + "_inv_pie_2020.png"
    plt.savefig(name)
    #plt.show()

