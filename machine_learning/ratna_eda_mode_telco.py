# -*- coding: utf-8 -*-
"""
Created on Mon Apr  7 13:32:57 2025

@author: ratsa
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


df = pd.read_csv('WA_Fn-UseC_-Telco-Customer-Churn.csv')

df.head(10) #first 10 rows
df.sample(20) #randomly picking

df.shape   #dimension
df.dtypes #datatypes


# Count how many times each data type is present in the dataset
pd.value_counts(df.dtypes)


df.info()
df.describe() #statistical summary

# TotalCharges is stored as an object; it needs conversion
df['TotalCharges'] =  pd.to_numeric(df['TotalCharges'], errors='coerce')

#or

# Replacing the blank values in the "TotalCharges" column with zeros and changing the data type to float.
df["TotalCharges"] = df["TotalCharges"].replace(" ", "0")



#printing range for TotalCharges 
range = df ['TotalCharges'].max() -df ['TotalCharges'].min()
print(f"Range: {range}")

#printing range for all numeric columns
range = df.max(numeric_only=True) -df.min(numeric_only=True)
print(f"Range: \n{range}")

df.columns

#display non numerical value
df.select_dtypes(exclude=["number","float","integer"]).head(10)

#checking for duplicate values, customer id is unique thats why it selected
df["customerID"].duplicated().sum()



# convert 'SeniorCitizen' to categorical
df['SeniorCitizen'] = df['SeniorCitizen'].map({1:'Yes', 0:'No'})



df.info()

#checking for null value
df.isnull().sum()

#drop columns with na
df.dropna(how='any', inplace=True)
df.info()





#Univariate analysis



df.hist(bins=10, figsize=(15, 10))
plt.show()


import plotly.express as px
import plotly.io as pio

# Set renderer to display inside Spyder
pio.renderers.default = "png"  # or use "png" for static plots

#boxplot
for column in df.select_dtypes(include=['float64', 'int64']).columns:
    fig = px.box(df, y=column, title=f'Box Plot of {column}')
    fig.show()
    

#how many customers churned out
plt.figure(figsize=(5,5))
abc = sns.countplot(x = 'Churn', data = df)
abc.bar_label(abc.containers[0])
plt.show()

# abc.containers[0]: In a countplot, the bars are organized into containers, and 
# abc.containers[0] refers to the first container that holds the bars for each category in 'Churn'.
# abc.bar_label(...): This method adds labels to the bars in the plot. 
# It takes the container of the bars and places the corresponding values (counts) on top of each bar.



# Creating pie chart of churned customers
plt.figure(figsize=(4,5))
gb = df.groupby('Churn').agg({'Churn':'count'})
plt.pie(gb['Churn'], labels = gb.index, autopct = "%1.2f%%")
plt.title("Percentage of Churned Customers", fontsize = 10)
plt.show()




# visualizing the distribution of each categorical feature with respect to churn
for idx, feature in enumerate(df.drop(columns=['Churn', 'TotalCharges', 'MonthlyCharges', 'tenure'])):
    plt.figure(idx, figsize=(8,5))
    sns.countplot(data=df, x=feature, hue='Churn')
    plt.title(f"Churn Distribution Across {feature}", fontweight='bold')
    plt.xticks(rotation=25)
    plt.legend(title='Churn Status')
    plt.show()
    print("\n\n")
    
#    Univariate Analysis:
#Customers who use streaming movies are less likely to churn.
#Month-to-month contracts have the highest churn.
#One-year and two-year contracts show very low churn rates.
#Churn is much higher for Fiber optic users compared to DSL users.
#Paperless billing users churn more compared to those who receive paper bills.
#Customers paying through electronic checks tend to churn the most.
#Automatic payments (bank transfer or credit card) users are more stable.
#New customers with less than a year of tenure churn the most.
#Customers with higher tenure rarely churn.



# visualizing the distribution of Monthly Charges for churned and non-churned customers
plt.figure(figsize=(8,6))
sns.kdeplot(df[df["Churn"]=='No']['MonthlyCharges'], shade=True, color='red', label="No Churn")
sns.kdeplot(df[df['Churn']=='Yes']['MonthlyCharges'], shade=True, color='blue', label="Churn")
plt.legend(loc='upper right')
plt.xlabel("Monthly Charges")
plt.title('Distribution of Monthly Charges Based on Churn Status')
plt.show()
#Insight: Customers with higher monthly charges tend to have a higher churn rate.


# visualizing the distribution of Total Charges for churned and non-churned customers
plt.figure(figsize=(8,6))
sns.kdeplot(df[df['Churn']=='No']['TotalCharges'], shade=True, color='red', label="No Churn")
sns.kdeplot(df[df['Churn']=='Yes']['TotalCharges'], shade=True, color='blue', label='Churn')
plt.legend(loc='upper right')
plt.xlabel('Total Charges')
plt.title("Distribution of Total Charges Based on Churn Status")
plt.show()


# visualizing the distribution of Tenure for churned and non-churned customers
plt.figure(figsize=(8,6))
sns.kdeplot(df[df['Churn']=='No']['tenure'], shade=True, color='red', label="No Churn")
sns.kdeplot(df[df['Churn']=='Yes']['tenure'], shade=True, color='blue', label='Churn')
plt.legend(loc='upper right', fontsize=20)  # Adjust the font size for the legend
plt.xlabel('Tenure', fontsize=20)  # Adjust the font size for the x-axis label
plt.ylabel('Tenure', fontsize=20) 
plt.title("Distribution of Total Charges Based on Churn Status", fontsize=20)  # Adjust the font size for the title
plt.show()




#bivariate analysis

# Relation berween Monthly Charges and Total Charges
sns.scatterplot(df, x='MonthlyCharges', y='TotalCharges')
plt.title("Total vs. Monthly Charges")
plt.show()


plt.figure(figsize=(8,5))
sns.boxplot(df,x="Churn", y='MonthlyCharges')
plt.title('Monthly Charges vs Churn')
plt.show()


#multivaruate analysis

# tenure + payment method vs churn
plt.figure(figsize=(10,6))
sns.boxplot(df, x='PaymentMethod', y='tenure', hue='Churn', showfliers=False)
plt.title('Tenure vs. Payment Method with respect to Churn', fontweight='bold')
plt.xlabel('Payment Method')
plt.ylabel('Tenure (Months)')
plt.xticks(rotation=45)
plt.legend(title='Churn', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.show()



# Get the correlation matrix and heatmap

# convert 'SeniorCitizen' to categorical


temp = df.copy()
temp['Churn'] = np.where(temp.Churn == 'Yes', 1, 0)
temp = temp.drop('tenure', axis=1)

#convert the categorical into dummy 0,1

temp1 = pd.get_dummies(temp, dtype=np.int32)

plt.figure(figsize=(12,12))
sns.heatmap(temp1.corr(), cmap='coolwarm', linewidths=0.5)
plt.title("Correlation Heatmap", fontweight='bold')
plt.show()

