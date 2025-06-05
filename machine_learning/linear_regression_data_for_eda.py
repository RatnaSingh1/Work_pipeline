# -*- coding: utf-8 -*-
"""
Created on Mon Apr  7 10:58:20 2025

@author: ratsa
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score

df = pd.read_csv('data_for_eda.csv')

#initialize the model
algorithm = LinearRegression() # here () different parameters can be use to optimize the model
model = algorithm.fit( X_train, ytrain)


# Define the features and the target variable
X = df.drop(['Income'], axis=1) #independent variables (predictors)
y = df['Income'] #dependent variables (target


# Split the data into training+validation and testing sets (80% training+validation, 20% testing)
X_train_val, X_test, y_train_val, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

X_train_val

y_train_val

y_test


# Split the training+validation set into training and validation sets (75% training, 25% validation of the 80%)
X_train, X_val, y_train, y_val = train_test_split(X_train_val, y_train_val, test_size=0.25, random_state=42)


#model Validation# Validate the model
y_val_pred = model.predict(X_val)
val_mse = mean_squared_error(y_val, y_val_pred) #mean squre error
accuracy = model.score(X_val, y_val) # accuracy (R-sqaured)

print(f"Validation Mean Squared Error: {val_mse}")

print(f"Validation Accuracy: {round(accuracy*100,2)}%")

#model testing

#  Make predictions on the test set
y_test_pred = model.predict(X_test)

# Evaluate the model

test_mse = mean_squared_error(y_test, y_test_pred) #mean squre error
accuracy = model.score(X_test, y_test) # accuracy (R-sqaured)



print(f"Test Mean Squared Error: {test_mse}")

print(f"Testing Accuracy: {round(accuracy*100,2)}%")



# Predict sample income
# Select a sample from the original dataframe
#n=1 one records selected randomly from df dataset stored in sample
sample = df.sample(n=1) # add random_state=42 parameter to get same results each time you run the code
print("Sample for prediction:")

# Preprocess the sample in the same way as the training data
sample_X = sample.drop(['Income'], axis=1)

sample_X
sample

# Select a sample from the original dataframe

# Convert the sample DataFrame to a NumPy array and reshape it
sample_X_array = sample_X.values.reshape(1, -1)  # Reshape to a 2D array with one row and the appropriate number of columns

# Predict the income for the sample
predicted_income = model.predict(sample_X_array)[0]  # Use the reshaped array for prediction
print(f"Predicted Income for the sample: {round(predicted_income,3)}")

# Display the actual income for the sample
actual_income = sample['Income'].values[0]
print(f"Actual Income for the sample: {round(actual_income,3)}")
print("Difference : ",round(predicted_income - actual_income,3))




#check for multiple samples
# Restore the built-in range function if necessary
if not callable(range):
    del range

# Create an empty list to store results
results = []

# Predict 5 different samples from the original df
for i in range(5):
    # Select a sample from the original dataframe
    sample = df.sample(n=1)

    # Preprocess the sample in the same way as the training data
    sample_X = sample.drop(['Income'], axis=1)
    sample_X_array = sample_X.values.reshape(1, -1)

    # Predict the income for the sample
    predicted_income = model.predict(sample_X_array)[0]

    # Store actual and predicted income in the results list
    results.append({
        "Sample": i+1,
        "Actual Income": sample['Income'].values[0],
        "Predicted Income": predicted_income,
        "Difference": sample['Income'].values[0] - predicted_income
    })
    
    # Create a DataFrame from the results
result_df = pd.DataFrame(results)
result_df
