# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 14:37:05 2025

@author: ratsa
"""

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from scipy import stats
from mpl_toolkits.mplot3d import Axes3D
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt # plotting
import os # accessing directory structure
import plotly.graph_objects as go

df = pd.read_csv("StudentsPerformance.csv")

#features

len(df)
df.info()
df.groupby(['gender']).size()
df.groupby(['race/ethnicity']).size()
df.groupby(['parental level of education']).size()
df.groupby(['test preparation course']).size()
df.describe()

df['math score'].describe()

# new column total score 
df['Total score']=df['math score']+df['reading score']+df['writing score']
df['Total score'].describe()



plt.figure(figsize=(18,8))
plt.subplot(1, 4, 1)
plt.title('MATH SCORES')
sns.barplot(y='math score',data=df,color='blue', ci='sd', linewidth=2)

# Adding mean annotations
mean_score = df['math score'].mean()
std_score = df['math score'].std()


#hue='gender': This splits the data by the 'gender' column, creating separate bars for each gender category.palette='seismic': This specifies the color palette to use for the bars. The 'seismic' palette is used here.
#\\{'hatch':'', 'alpha':0.6, 'linewidth':2}*:hatch=''*: This adds a hatch pattern (a star) to the bars.alpha=0.6: This sets the transparency of the bars to 60% (0.6).
#linewidth=2: This sets the width of the lines outlining the bars to 2.

plt.figure(figsize=(14,8))
plt.subplot(1, 3, 1)
sns.barplot(x='test preparation course',y='math score',data=df,hue='gender',palette='seismic',**{'hatch':'*','alpha':0.6,'linewidth':2})
plt.title('MATH SCORES')

#histo
f,(ax_box, ax_hist) = plt.subplots(2,figsize=(16,8),sharex=True,gridspec_kw={"height_ratios": (.50, .85)})
sns.boxplot(x='math score',data=df,color='#ff1a1a',ax=ax_box)
sns.histplot(data = df,x = 'math score',bins=20,edgecolor='black',color='#ff1a1a',kde=True,ax=ax_hist)
plt.title('Math score data distribution',color='black',size=25)
plt.show()

#gender pie chart
df = pd.read_csv("StudentsPerformance.csv")
gender_counts = df['gender'].value_counts()
fig = go.Figure([go.Pie(labels=gender_counts.index, values=gender_counts, opacity=0.9)])
fig.update_traces(textinfo='percent+label', marker=dict(line=dict(color='#000000', width=2)))
fig.update_layout(title_text='Distribution of the Gender', title_x=0.5, title_font=dict(size=20))
fig.show()


import pandas as pd
import plotly.graph_objects as go  # Ensure this is imported

# Load the dataset
df = pd.read_csv("StudentsPerformance.csv")

# Get the gender counts
# Get the gender counts
gender_counts = df['gender'].value_counts()
# Create subplots
f, ax = plt.subplots(1, 2, figsize=(20, 10))
# Pie chart
ax[0].pie(gender_counts, labels=gender_counts.index, autopct="%1.1f%%", startangle=90)
ax[0].set_title('Gender Distribution (Pie Chart)')
# Count plot
sns.countplot(x='gender', data=df, ax=ax[1])
ax[1].set_title('Gender Distribution (Count Plot)')
# Display the plots
plt.tight_layout()
plt.show()






plt.subplot(1, 3, 3)
sns.barplot(x='test preparation course',y='writing score',data=df,hue='gender',palette='seismic',**{'hatch':'x','linewidth':2})
plt.title('WRITING SCORES')
plt.show()

#or

plt.rcParams['figure.facecolor'] = "#e6ecff"
plt.rcParams['axes.facecolor'] = "#e6ecff"
plt.figure(figsize=(14,8))
plt.subplot(1, 3, 1)
sns.barplot(x='test preparation course',y='math score',data=df,hue='gender',palette='seismic',**{'hatch':'*','alpha':0.6,'linewidth':2})
plt.title('MATH SCORES')
plt.subplot(1, 3, 2)
sns.barplot(x='test preparation course',y='reading score',data=df,hue='gender',palette='seismic',**{'hatch':'.','alpha':0.8,'linewidth':2})
plt.title('READING SCORES')
plt.subplot(1, 3, 3)
sns.barplot(x='test preparation course',y='writing score',data=df,hue='gender',palette='seismic',**{'hatch':'x','linewidth':2})
plt.title('WRITING SCORES')
plt.show()




#hypothesis testing

import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt

#Sample weekly sales data (before and after product launch)

import pandas as pd
from scipy import stats

# Weekly sales data
weekly_sales = {
    'pre_launch': [1000, 850, 920, 1100, 880],
    'post_launch': [1200, 1350, 1180, 1500, 1420]}

# Load data into a pandas DataFrame
df = pd.DataFrame(weekly_sales)

# Calculate average sales before and after launch
pre_launch_mean = df['pre_launch'].mean()
post_launch_mean = df['post_launch'].mean()

# Perform t-test
tstatistic, pvalue = stats.ttest_ind(df['pre_launch'], df['post_launch'])

# Print results
print(f"Average weekly sales before launch: {pre_launch_mean:.2f}")
print(f"Average weekly sales after launch: {post_launch_mean:.2f}")

if pvalue < 0.05:
    print(f"There is a statistical difference in sales (p-value:", pvalue,")")
else:
    print(f"There is no statistical difference in sales (p-value:", pvalue,")")

#or

if pvalue < 0.05:
    print(f"There is a statistical difference in sales (p-value: {pvalue:.3f})")
else:
    print(f"There is no statistical difference in sales (p-value: {pvalue:.3f})")
    
    

# Create a bar plot
means = [pre_launch_mean, post_launch_mean]
labels = ['Pre-Launch', 'Post-Launch']

plt.bar(labels, means, color=['blue', 'orange'])
plt.ylabel('Average Sales')
plt.title('Weekly Sales Before and After Launch')

# Add significance annotation
x1, x2 = 0, 1  # x-coordinates of the two groups
y, h = max(means) + 50, 50  # y-coordinate for the line and the height
plt.plot([x1, x2], [y, y], color='black', linewidth=1)  # Line connecting the groups
plt.text((x1 + x2) / 2, y + h, f'p = {pvalue:.3f}', ha='center', va='bottom', color='black')

# Display the plot
plt.show()



#another way of plotting

import pandas as pd 
import plotly.graph_objects as go

# Data

data = {
    'pre_launch sales': [1000, 850, 920, 1100, 880],
    'post_launch sales': [1200, 1350, 1180, 1500, 1420]}

#create data frame

df = pd.DataFrame(data, index=['Product 1', 'Product 2', 'Product 3', 'Product 4', 'Product 5' ])

# create  bar plot using plotly

fig = go.Figure()

fig.add_trace(go.Bar(
    x=df.index,
    y=df['pre_launch sales'],
    name='Pre-Launch Sales',
    marker_color='blue'   
    ))


fig.add_trace(go.Bar(
    x=df.index,
    y=df['post_launch sales'],
    name='Post-Launch Sales',
    marker_color='seagreen'   
    ))

#add labels , title, and legend

fig.update_layout(
    title='Comparision of Pre-Launch and Post-Launch Sales',
    xaxis=dict(title='Products'),
    yaxis=dict(title='Sales'),
    barmode = 'group',
    width=1000,
    height=1000 )
fig.show()
    
