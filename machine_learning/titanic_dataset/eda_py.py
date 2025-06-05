# -*- coding: utf-8 -*-
"""
Created on Tue Apr  1 15:44:08 2025

@author: ratsa
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


df = pd.read_csv('train.csv')

df.head(10) #first 10 rows
df.sample(20) #randomly picking

df.shape   #dimension
df.dtypes #datatypes


# Count how many times each data type is present in the dataset
pd.value_counts(df.dtypes)

df.info()
df.describe() #statistical summary

df = df.drop('Cabin', axis=1)
df = df.drop('Name', axis=1)

df.columns

#display non numerical value
df.select_dtypes(exclude=["number","float","integer"]).head(10)


import matplotlib.pyplot as plt

df.hist(bins=20, figsize=(15, 10))
plt.show()




import plotly.express as px
import plotly.io as pio

# Set renderer to display inside Spyder
pio.renderers.default = "png"  # or use "png" for static plots

#boxplot
for column in df.select_dtypes(include=['float64', 'int64']).columns:
    fig = px.box(df, y=column, title=f'Box Plot of {column}')
    fig.show()
    
    
#bivariate analysis
    
    df.columns
    
    import plotly.express as px

list = ['Fare','Survived']

fig = px.scatter(df, x=list[0], y=list[1], title=f'Scatter Plot of {list[0]} vs {list[1]}')
fig.show()


list =['Fare','Survived']

fig = px.scatter(df, x=list[0], y=list[1], title=f'Scatter Plot of {list[0]} vs {list[1]}')
fig.show()


columns = ['Age']

df2 = df.dropna(subset=columns)  # Drops rows with NaN in 'Age' or 'Income'
fig = px.scatter(df2, x=columns[0], y=columns[1], title=f'Scatter Plot of {columns[0]} vs {columns[1]}')
fig.show()





import plotly.express as px

for column in df.select_dtypes(include=['object']).columns:
    fig = px.bar(df, x=column, title=f'Bar Chart of {column}')
    fig.show()
    
    
import matplotlib.pyplot as plt
import seaborn as sns

# Get the correlation matrix and heatmap


numeric_df = df.select_dtypes(include=['float64', 'int64'])
corr = numeric_df.corr()
corr


# Set the figure size
plt.figure(figsize=(10, 8))

# Generate the heatmap
sns.heatmap(corr, annot=True, cmap='coolwarm')

# Display the heatmap
plt.show()



#multivariate analysis pairplot

# plot a pairplot for df

import seaborn as sns

# Create a pairplot using Seaborn
sns.pairplot(df)
plt.show()


df_filtered = df[df['Fare'] > 500]
df_filtered



# Swarm Plot: Fare vs Survival
sns.swarmplot(x="Survived", y="Fare",hue="Pclass",data=df, palette="coolwarm")
plt.title("Fare vs Survival")
plt.xlabel('Survival')
plt.show()

#pie plot survived vs not survived
df["Survived"].value_counts().plot.pie(autopct="%1.1f%%", labels=["Not Survived", "Survived"], colors=["pink", "skyblue"],shadow=True)
plt.title("Overall Survival Percentage")
plt.ylabel("")  # Remove y-label for aesthetics
plt.show()


# Bar Plot: Survival Rate by Class
plt.figure(figsize=(8, 5))
ax=sns.barplot(x=df["Pclass"], y=df["Survived"]*100,hue=df['Sex'],ci=None,palette="dark:#5A9_r")
# Add value labels on top of bars
for container in ax.containers:
    ax.bar_label(container, fmt="%.1f%%", fontsize=9)  # Format as percentage
plt.xlabel("Passenger Class",fontsize=12)
plt.ylabel("Survival Rate %",fontsize=12)
plt.title("Survival Rate %age by Passenger Class",fontsize=14)
plt.xticks(fontsize=12)
plt.show()


#making a age group

df['Age'].value_counts().sort_index()

def categorize_age(age):
    if age <= 4:
        return "Baby"
    elif age <= 12 and age>=5:
        return "Child"
    elif age>=13 and age <= 19:
        return "Teen"
    elif age>=20 and age <= 39:
        return "Adult"
    elif age >=40 and age <= 59:
        return "Middle Age Adult"
    else:
        return "Senior Adult"

df["Age_group"] = df["Age"].apply(categorize_age)
df["Age_group"].value_counts()
