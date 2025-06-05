# -*- coding: utf-8 -*-
"""
Created on Mon Mar 31 13:20:59 2025

@author: ratsa
"""

import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report
import warnings


df = pd.read_csv('hypothyroid.csv')
df.head(10)


# Display basic information about the dataset
print(df.info())


# Display summary statistics
df.describe()
df.dtypes


df[['age', 'TSH', 'T3', 'TT4', 'T4U', 'FTI']]

# Check if "?" exist in the df. if exist, show them

df[df.applymap(lambda x: x == '?').any(axis=1)] # Use axis=1 to check if any value in each row is True

# Replace '?' with NaN
df.replace('?', np.nan, inplace=True)
df[df.isnull().any(axis=1)]
df.info()




# Convert columns to appropriate data types
numeric_columns = ['age', 'TSH', 'T3', 'TT4', 'T4U', 'FTI']
for col in numeric_columns:
    df[col] = pd.to_numeric(df[col], errors='coerce')
    
    df[numeric_columns]

df.dtypes





#preprocess the data 
#dealing outliers


import plotly.express as px

for col in numeric_columns:
    fig = px.box(df, y = col, title=f"Boxplot of {col}")
    fig.show()


#removing outliers using quartile


print(f"DataFrame dimensions BEFORE removing outliers: {df.shape}")

# Select only numerical columns for outlier removal
numerical_df = df.select_dtypes(include=['float64', 'int64'])

Q1 = numerical_df.quantile(0.25)
Q3 = numerical_df.quantile(0.75)
IQR = Q3 - Q1

# Apply outlier removal to the numerical DataFrame
numerical_df = numerical_df[~((numerical_df < (Q1 - 1.5 * IQR)) |(numerical_df > (Q3 + 1.5 * IQR))).any(axis=1)]

# If you need to update the original df, you can reassign the numerical columns
# df[numerical_df.columns] = numerical_df

print(f"DataFrame dimensions AFTER removing outliers: {numerical_df.shape}")

numerical_df

df[numerical_df.columns] = numerical_df
# plot a boxplot for each numeric feature
for col in numeric_columns:
    fig = px.box(df, y=col, title=f"Boxplot of {col}")
    fig.show()

#missing value

df.isnull().sum()

# Import the necessary library


import plotly.express as px
import plotly.io as pio

# Set renderer to display inside Spyder
pio.renderers.default = "png"  # or use "png" for static plots

# Plot the distribution of each numeric feature
for col in numeric_columns:
    fig = px.histogram(df, x=col, nbins=10, title=f"Distribution of {col}")
    fig.show()


#replacing values with mean to get normal distribution

#checking max and minage and then mean
df['age'].max()

df['age'].min()
df['age'].mean()

#imputing by mean which have nulll value

from sklearn.impute import SimpleImputer

imputer = SimpleImputer(strategy='mean')
df[['age','T3', 'TT4', 'T4U', 'FTI']] = imputer.fit_transform(df[['age','T3', 'TT4', 'T4U', 'FTI']])  # Replace with appropriate columns

#NB: age has only 1 extreme value (see the boxplot)

df.isnull().sum()

#replacing the missing value with median in TSH
#fill the missing values in TSH column with the median

df['TSH'].fillna(df['TSH'].median(), inplace=True)

df.isnull().sum()



# Fill missing values for categorical columns with the mode
categorical_columns = df.select_dtypes(include=['object']).columns
for col in categorical_columns:
    df[col].fillna(df[col].mode()[0], inplace=True)

df.isnull().sum()


# Drop columns with all missing values (if any)
df.dropna(axis=1, how='all', inplace=True)


df.isnull().sum()


#Feature engineering

df['on thyroxine'].nunique()

df.columns

#count the unique features in each column and put the results in a dataframe. also show the particular unique features present

unique_features = {}
for col in df.columns:
  unique_features[col] = {
      "count": df[col].nunique(),
      "unique_values": df[col].unique()
  }

unique_features_df = pd.DataFrame(unique_features).T
unique_features_df


df.columns



# Convert binary categorical columns to 0 and 1
binary_columns = [
    "sex", "on thyroxine", "query on thyroxine", "on antithyroid medication", "sick",
    "pregnant", "thyroid surgery", "I131 treatment", "query hypothyroid", "query hyperthyroid",
    "lithium", "goitre", "tumor", "hypopituitary", "psych", "TSH measured", "T3 measured",
    "TT4 measured", "T4U measured", "FTI measured"
]

#checking after transformation
for col in binary_columns:
    df[col] = df[col].map({'t': 1, 'f': 0, 'M': 1, 'F': 0})

unique_features = {}
for col in df.columns:
  unique_features[col] = {
      "count": df[col].nunique(),
      "unique_values": df[col].unique()
  }
  
  unique_features_df = pd.DataFrame(unique_features).T
  unique_features_df


# convert TBG measured class from "f" to 1

df['TBG measured'].replace('f', 1, inplace=True)


# @title
#count the unique features in each column and put the results in a dataframe. also show the particular unique features present

unique_features = {}
for col in df.columns:
  unique_features[col] = {
      "count": df[col].nunique(),
      "unique_values": df[col].unique()
  }

unique_features_df = pd.DataFrame(unique_features).T
unique_features_df




#convert the classes in  "referral source" feature to numeric classes

referral_source_mapping = {
    'SVHC': 0,
    'other': 1,
    'SVI': 2,
    'STMW': 3,
    'SVHD': 4,
}

df['referral source'] = df['referral source'].map(referral_source_mapping)

# @title
#count the unique features in each column and put the results in a dataframe. also show the particular unique features present

unique_features = {}
for col in df.columns:
  unique_features[col] = {
      "count": df[col].nunique(),
      "unique_values": df[col].unique()
  }

unique_features_df = pd.DataFrame(unique_features).T
unique_features_df




# Convert the target variable to binary
df['binaryClass'] = df['binaryClass'].map({'P': 1, 'N': 0})

# @title
#count the unique features in each column and put the results in a dataframe. also show the particular unique features present

unique_features = {}
for col in df.columns:
  unique_features[col] = {
      "count": df[col].nunique(),
      "unique_values": df[col].unique()
  }

unique_features_df = pd.DataFrame(unique_features).T
unique_features_df
df
df.types

#till here everything changed to numeric data thus ready for machine learning model

#splitting the data
#building the model, drop the binary class( to predict)

# Define the features and the target variable
X = df.drop(columns=['binaryClass'])    #features
y = df['binaryClass']                      #target
df
y.value_counts()

#splitting
from sklearn.model_selection import train_test_split



# Split the data into training+validation and testing sets (80% training+validation, 20% testing)
X_train_val, X_test, y_train_val, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Split the training+validation set into training and validation sets (75% training, 25% validation of the 80%)
X_train, X_val, y_train, y_val = train_test_split(X_train_val, y_train_val, test_size=0.25, random_state=42)


#standerization after splitting to prevent data leakage, only standarizied the testing data

# @title
from sklearn.preprocessing import StandardScaler

# Initialize the scaler
scaler = StandardScaler()

# Fit the scaler on the training data and transform the training and validation data
X_train_scaled = scaler.fit_transform(X_train)
X_val_scaled = scaler.transform(X_val)
X_test_scaled = scaler.transform(X_test)





#model training

# Train the model
algorithm= LogisticRegression(max_iter=1000, random_state=42)
model = algorithm.fit(X_train_scaled, y_train)


# Validate the model
y_val_pred = model.predict(X_val_scaled)
val_accuracy = accuracy_score(y_val, y_val_pred)
# val_conf_matrix = confusion_matrix(y_val, y_val_pred)
# val_class_report = classification_report(y_val, y_val_pred)
print(f"Validation Accuracy: {round(val_accuracy*100,2)}%")

df.binaryClass.value_counts() # in output we have more 1 than 0 , that is called data imbalance


# Evaluate the model on the test set
y_test_pred = model.predict(X_test_scaled)
test_accuracy = accuracy_score(y_test, y_test_pred)
# test_conf_matrix = confusion_matrix(y_test, y_test_pred)
# test_class_report = classification_report(y_test, y_test_pred)
print(f"Test Accuracy: {round(test_accuracy*100,2)}%")


#model prediction
#select random 10 samples from index

#
sample = X_test.sample(n=10)
sample

index = sample.index.to_list()
index

 sample_indices = index
 X_test.loc[sample_indices]
 
 
 # Make sample predictions, store in sample data, and use scaler to transform sample data
sample_indices = index  # Indices of samples to predict
sample_data = X_test.loc[sample_indices]
sample_data_scaled = scaler.transform(sample_data)             #scaling
sample_predictions = model.predict(sample_data_scaled)          #prediction

# Display sample predictions
print("\nSample Predictions:")
for i, idx in enumerate(sample_data.index): # Iterate over the index of the sample_data DataFrame
    print(f"Sample Index: {idx}")
    print(f"Features: {sample_data.iloc[i].to_dict()}") # Use iloc to index relative to sample_data
    print(f"Actual Class: {y_test.loc[idx]}") # Use loc with the original index to access y_test
    print(f"Predicted Class: {sample_predictions[i]}")
    print()    #check the prediction here
    
    
    #
# put the result in a dataframe in the "Sample" dataframe where we selected the samples at first

predictions_df = pd.DataFrame({
    "Sample Index": sample_indices,
    "Features": sample_data.to_dict(orient="records"),
    "Actual Class": y_test.loc[sample_indices],
    "Predicted Class": sample_predictions
})

predictions_df
    
    # feature importance
    
#model coeffient higherthe coefficient higher the contribution of that feature, this 
#give coeffiecient of features we have used    
model.coef_[0]



#
import plotly.express as px

# Get feature importances
importance = np.abs(model.coef_[0])  #calculating absolute value of model_coef(due to -ve,+ve values)
feature_names = X.columns  # extracting all columns except binary class

# Create a DataFrame for visualization from selected coloumns and coef
feature_importance_df = pd.DataFrame({'Feature': feature_names, 'Importance': importance})

# Sort by importance
feature_importance_df = feature_importance_df.sort_values(by='Importance', ascending=True)

# Plot the feature importances
fig = px.bar(feature_importance_df, x='Importance', y='Feature', orientation='h', color_discrete_sequence=['#008080'], title='Feature Importances')
fig.update_layout(
    plot_bgcolor='rgb(17, 17, 17)',
    paper_bgcolor='rgb(17, 17, 17)',
    font_color='white'
)
fig.show()






# Select top N features (e.g., top 10 features) for modelbuilding
top_features = feature_importance_df.head(10)['Feature'].values

# Create a new dataset with only the top features
X_top = X[top_features]

# Split the data into training+validation and testing sets
X_train_val_top, X_test_top, y_train_val, y_test = train_test_split(X_top, y, test_size=0.2, random_state=42)
X_train_top, X_val_top, y_train, y_val = train_test_split(X_train_val_top, y_train_val, test_size=0.25, random_state=42)

# Standardize the data
X_train_top_scaled = scaler.fit_transform(X_train_top)
X_val_top_scaled = scaler.transform(X_val_top)
X_test_top_scaled = scaler.transform(X_test_top)

# Train the model with the top features
model_top = LogisticRegression(penalty='l2', C=1.0, max_iter=1000, random_state=42)
model_top.fit(X_train_top_scaled, y_train)

# Validate the model
y_val_pred_top = model_top.predict(X_val_top_scaled)
val_accuracy_top = accuracy_score(y_val, y_val_pred_top)
val_conf_matrix_top = confusion_matrix(y_val, y_val_pred_top)
val_class_report_top = classification_report(y_val, y_val_pred_top)
print(f"Validation Accuracy with Top Features: {round(val_accuracy_top*100,2)}%")
# print(f"Validation Confusion Matrix with Top Features:\n{val_conf_matrix_top}")
# print(f"Validation Classification Report with Top Features:\n{val_class_report_top}")

# Evaluate the model on the test set
y_test_pred_top = model_top.predict(X_test_top_scaled)
test_accuracy_top = accuracy_score(y_test, y_test_pred_top)
test_conf_matrix_top = confusion_matrix(y_test, y_test_pred_top)
test_class_report_top = classification_report(y_test, y_test_pred_top)
print(f"Test Accuracy with Top Features: {round(test_accuracy_top*100,2)}%")
# print(f"Test Confusion Matrix with Top Features:\n{test_conf_matrix_top}")
# print(f"Test Classification Report with Top Features:\n{test_class_report_top}")

# Make sample predictions
sample_indices = [0, 10, 20]  # Indices of samples to predict
sample_data_top = X_test_top.iloc[sample_indices]
sample_data_top_scaled = scaler.transform(sample_data_top)
sample_predictions_top = model_top.predict(sample_data_top_scaled)








# # Display sample predictions
# print("\nSample Predictions with Top Features:")
# for i, idx in enumerate(sample_indices):
#     print(f"Sample Index: {idx}")
#     print(f"Features: {X_test_top.iloc[idx].to_dict()}")
#     print(f"Actual Class: {y_test.iloc[idx]}")
#     print(f"Predicted Class: {sample_predictions_top[i]}")
#     print()
