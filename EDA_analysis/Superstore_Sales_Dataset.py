# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 10:39:33 2025

@author: ratsa
"""

#read file

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.tsa.arima_model import ARIMA
from statsmodels.tsa.seasonal import seasonal_decompose

data = pd.read_csv('train.csv')

data.head()
data.shape
data.columns
data.info()
data.describe()
data['Sales'].describe()


#Setting up row ID as the index
data.set_index("Row ID", inplace=True)

#handling missing values
data.isnull()

data.isnull().sum()  # Check for missing values

#data = data[data.isnull().any(axis=1)]

# Fill missing values with a specific value
data = data.fillna(0)
#removing duplicates

data = data.drop_duplicates()

data = data.dropna()


#save data
data.to_csv('cleaned_data.csv', index=False)  # Save to CSV file

#grouping

data.groupby(['Segment'])['Sales'].mean()



# Add a seperate column year to the table
data["Order Date"]=pd.to_datetime(data["Order Date"], format='%d/%m/%Y')
data["Year"]=data["Order Date"].dt.year

data.head()


# sorting
data1= data.sort_values('Year')


# Group sales by year and mean sales
                 
  data.groupby(['Year'])['Sales'].mean()               
                 
                

# Group sales by state and sum them
sales_by_state = data.groupby('State')['Sales'].sum().sort_values()

# Create a bar chart for sales by state
plt.figure(figsize=(12,8))
plt.barh(sales_by_state.index, sales_by_state.values, color='skyblue')
plt.title('Total Sales by State')
plt.xlabel('Total Sales Amount')
plt.ylabel('State')
plt.grid(True, axis='x')
plt.tight_layout()
plt.show()

                 




# Convert the 'Order Date' column to datetime format with the specified format
data['Order Date'] = pd.to_datetime(data['Order Date'], format='%m/%d/%Y', errors='coerce')

# Verify the changes
print(data.head())



# Standardize categorical variables by converting to lowercase
data['Category'] = data['Category'].str.lower()

# Alternatively, you can convert to uppercase if needed
#data['Category'] = data['Category'].str.upper()

# Remove leading and trailing whitespaces
data['Category'] = data['Category'].str.strip()




#to get summary of sales
# Calculate summary statistics
summary = data['Sales'].agg(['sum', 'mean', 'count'])

# Rename the summary index for clarity
summary.index = ['Total Sales', 'Average Sales', 'Order Count']

print(summary)
# sales by category
# Group by 'Category' and calculate the total sales for each category
sales_by_category = data.groupby('Category')['Sales'].mean()

# Display the results
print(sales_by_category)

plt.figure(figsize=(12,8))
plt.barh(sales_by_category.index, sales_by_category.values, color='skyblue')
plt.show()

#sales by region and top 5 sales state to do

top_5_states = data.groupby('Region')['Sales'].sum().nlargest(5).reset_index()

print("Top 5 Sales States:")
print(top_5_states)






#sales trend anlysis

# Convert the 'Order Date' column to datetime format
data['Order Date'] = pd.to_datetime(data['Order Date'], format='%m/%d/%Y')

# Set the 'Order Date' column as the index
data.set_index('Order Date', inplace=True)

# Resample the data by month and calculate the total sales for each month
monthly_sales = data['Sales'].resample('M').sum()

# Plot the sales trend
plt.figure(figsize=(10, 6))
monthly_sales.plot()
plt.title('Monthly Sales Trend')
plt.xlabel('Date')
plt.ylabel('Total Sales')
plt.show()


#time series analysis


# Load your DataFrame
data = pd.read_csv('train.csv') 

# Convert the 'Order Date' column to datetime format
data['Order Date'] = pd.to_datetime(data['Order Date'], format='%d/%m/%Y')

# Set the 'Order Date' column as the index
data.set_index('Order Date', inplace=True)

# Resample the data by month and calculate the total sales for each month
monthly_sales = data['Sales'].resample('ME').sum()

# Decompose the time series into trend, seasonal, and residual components
decomposition = seasonal_decompose(monthly_sales, model='additive')
decomposition.plot()
plt.show()

# Fit an ARIMA model
model = ARIMA(monthly_sales, order=(5, 1, 0))
model_fit = model.fit(disp=0)

# Forecast future sales
forecast = model_fit.forecast(steps=12)[0]

# Plot the sales trend and forecast
plt.figure(figsize=(10, 6))
plt.plot(monthly_sales, label='Actual Sales')
plt.plot(pd.date_range(start=monthly_sales.index[-1], periods=12, freq='M'), forecast, label='Forecasted Sales', color='red')
plt.title('Sales Trend and Forecast')
plt.xlabel('Date')
plt.ylabel('Total Sales')

plt.legend()
plt.show()


# Add a seperate column year to the table
data["Order Date"]=pd.to_datetime(data["Order Date"], format='%d/%m/%Y')
data["Year"]=data["Order Date"].dt.year

data.head()

print(data["Year"].max())
print(data["Year"].min())

#Grouping the sales with year and creating it in a new dataframe
data_SY=data.groupby("Year")["Sales"].mean().reset_index()
data_SY
sns.barplot(data=data_SY,x="Year",y="Sales", color="violet")
plt.xlabel("Year")
plt.ylabel("Yearly sales")
plt.title("Sales of  Per year")


# creating a dataset showing ftype of category highest sale
#Now creating a dataframe that shows the type of furniture are selling per year
data_TYF=data.groupby(["Category"])["Sales"].mean().reset_index()
data_TYF

sns.barplot(data=data_TYF, x="Category", y="Sales", color="red")
#plt.figure(figsize=(10,8))
plt.xticks(rotation=45)
plt.show()


#Now creating a dataframe that shows the type of furniture are selling per year
data_TYF=data.groupby(["Sub-Category"])["Sales"].mean().reset_index()
data_TYF

sns.barplot(data=data_TYF, x="Sub-Category", y="Sales", color="red")
#plt.figure(figsize=(10,8))
plt.xticks(rotation=45)
plt.show()




#Now we can do a pie chart to see which country has used which type of shipmode
data_SHM=data.groupby(["Country","Ship Mode"]).size().reset_index(name="Count")
data_SHM.set_index("Country", inplace=True)
data_SHM


#Now make a Pie chart to see the variation
data_SHM["Count"].plot(kind="pie",startangle=90, figsize=(5,6),autopct='%1.1f%%', legend=True, ylabel='', labels=data_SHM["Ship Mode"])
plt.title("Different categories of Shipping mode used by US")
plt.show()




# Group by Category and Year to get total sales for each combination
#Pivots the DataFrame to reshape it, making 'Year' the index and 'Category' the columns.
sales_by_category_year = data.groupby(['Category', 'Year'])['Sales'].sum().reset_index()
# Pivot the data to reshape it
pivot_table = sales_by_category_year.pivot(index='Year', columns='Category', values='Sales')

# Plot the data
pivot_table.plot(kind='bar', figsize=(10, 6))

plt.xlabel('Year')
plt.ylabel('Sales')
plt.title('Sales by Category and Year')
plt.legend(title='Category')


#sales by region
sales_by_category_region = data.groupby(['Region', 'Category'])['Sales'].sum().reset_index()
pivot_table1 = sales_by_category_region.pivot(index='Region', columns='Category', values='Sales')
pivot_table1.plot(kind='bar', figsize=(10, 6))



#Which city gave the highesct sale
data_CS=data.groupby("Category")["City"].size().reset_index()
data_CS
data_CS["City"].plot(kind="pie", startangle=90, figsize=(5,6),autopct='%1.1f%%', legend=True, ylabel='', labels=data_CS["Category"])
plt.title("Different categroeis according to different cities")
plt.show()




#We plot a histogram to see how many numbers of sub categories are in each category
data_grouped = data.groupby(["Category", "Sub-Category"]).size().unstack()
# Plot as a bar chart
data_grouped.plot(kind="bar", figsize=(8,6), stacked=True)
# Add title
plt.title("Sub-Category Frequency per Category")
plt.ylabel("Count")
plt.show()

# value counts

data['Ship Mode'].value_counts()
data['Ship Mode'].value_counts().plot.pie()

data["Segment"].value_counts().plot.pie()


# top 10 stats
top_stats = data.groupby(['State']).sum().sort_values('Sales', ascending=False).head(10)
top_stats = top_stats[['Sales']].round(2)
top_stats.reset_index(inplace=True)
top_stats




# What is the distribution of customers across different regions?

#firstly, identifying the unique customers
customers_per_region = data.groupby('Region')['Customer Name'].nunique().reset_index()

#now, visualizing the data
sns.barplot(data=customers_per_region, x='Region', y='Customer Name')
plt.ylabel('Unique Customers')
plt.xlabel('Regions')
plt.show()




# Which categories and sub-categories are the most profitable

#firstly, calculating the total and average profit for each category and sub-category
profit = data.groupby(['Category', 'Sub-Category'], as_index=False).agg({
    'Sales': ['sum', 'mean']
})

#renaming the columns name
profit.columns = ['category', 'sub_category', 'total_profit', 'average_profit']

#sorting by total profit
result = profit.sort_values(by='total_profit', ascending=False).head(10)

#visualizing the data
sns.barplot(data=result, x='total_profit', y='sub_category', hue='category')
plt.show()





