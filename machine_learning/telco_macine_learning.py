# -*- coding: utf-8 -*-
"""
Created on Thu Apr 10 15:47:40 2025

@author: ratsa
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Apr  7 10:58:20 2025

@author: ratsa
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import LabelEncoder

from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score

df = pd.read_csv('train.csv')
df.info()


#drop missing value
df=df.drop(['Cabin'], axis=1) 
df=df.drop(['Name'], axis=1) 
df=df.drop(['PassengerID'], axis=1) 
df['Age'].mean()
# replace the missing value in the Age feature with the mean

df['Age'].fillna(df['Age'].mean(), inplace=True)

#drop columns with na
df.dropna(how='any', inplace=True)
df.info()


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






# label encoding for binary columns
binary_cols = ['Sex', 'Ticket', 'Survived', 'SibSp', 'Fare' ]
le = LabelEncoder()
for col in binary_cols:
    df[col] = le.fit_transform(df[col])


# encoding remaining multiclass features using one hot encoding
multiclass_cols = ['Pclass', 'Embarked', ]

df = pd.get_dummies(df, columns=multiclass_cols, drop_first=True, dtype=np.int32)

print(df.isna().sum()) # testing column with na



# or use imputer
#from sklearn.impute import SimpleImputer
#imputer = SimpleImputer(strategy='mean')
#df['TotalCharges'] = imputer.fit_transform(df[['TotalCharges']])



df.head()

# Normalize the data
#scaler = MinMaxScaler()
#df = pd.DataFrame(scaler.fit_transform(df), columns=df.columns)







# Define the features and the target variable
X = df.drop(['Survived'], axis=1) #independent variables (predictors) 
y = df['Survived'] #dependent variables (target


from collections import Counter
Counter(y)

## handling imbalanced data
from imblearn.combine import SMOTEENN
sm = SMOTEENN(random_state=0)
X_resampled, y_resampled = sm.fit_resample(X,y)

Counter(y_resampled)

#print(X_train.isna().sum()) # testing column with na

#from sklearn.impute import SimpleImputer

#imputer = SimpleImputer(strategy='mean')  # Use 'median' if it's more appropriate
#df['TotalCharges'] = imputer.fit_transform(df[['TotalCharges']])




# Split the data into training+validation and testing sets (80% training+validation, 20% testing)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

X_train

y_train

y_test



from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.linear_model import LogisticRegression,LogisticRegressionCV,SGDClassifier
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier,AdaBoostClassifier,BaggingClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier


# Define models
models = {
    'Logistic Regression': LogisticRegression(),
        'Logistic Regression CV': LogisticRegressionCV(),
    'SGD': SGDClassifier(),
    'Random Forest': RandomForestClassifier(),
    'Gradient Boosting': GradientBoostingClassifier(),
    'AdaBoost': AdaBoostClassifier(),
    'Bagging': BaggingClassifier(),
    'Decision Tree': DecisionTreeClassifier(),
    'Support Vector Machine': SVC(),
    'K-Nearest Neighbors': KNeighborsClassifier()
}

# Train and evaluate models
def evaluate_models(X_train, X_test, y_train, y_test):
    results = []
    for name, model in models.items():
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)
        acc = accuracy_score(y_test, y_pred)
        results.append((name, acc))
    
    # Sort models by accuracy
    results.sort(key=lambda x: x[1], reverse=True)
    return results


results = evaluate_models(X_train,X_test,y_train,y_test)
    
print("Model Performance:")
for name, acc in results:
    print(f"{name}: {acc:.6f}")




#initialize the model
algorithm = LinearRegression() # here () different parameters can be use to optimize the model

model = algorithm.fit( X_train, y_train)






#model Validation# Validate the model
y_val_pred = model.predict(X_test)
val_mse = mean_squared_error(X_test, y_val_pred) #mean squre error



accuracy = model.score(X_train, y_val_pred) # accuracy (R-sqaured)

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
sample_X = sample.drop(['Churn'], axis=1)

sample_X
sample

# Select a sample from the original dataframe

# Convert the sample DataFrame to a NumPy array and reshape it
sample_X_array = sample_X.values.reshape(1, -1)  # Reshape to a 2D array with one row and the appropriate number of columns

# Predict the income for the sample
predicted_Churn = model.predict(sample_X_array)[0]  # Use the reshaped array for prediction
print(f"Predicted Churn for the sample: {round(predicted_Churn,3)}")

# Display the actual income for the sample
actual_income = sample['Churn'].values[0]
print(f"Actual Income for the sample: {round(actual_income,3)}")
print("Difference : ",round(predicted_income - actual_income,3))




#check for multiple samples
# Restore the built-in range function if necessary
if not callable(range):
    del range

# Create an empty list to store results
results = []

# Predict 5 different samples from the original df
for i in range(1000):
    # Select a sample from the original dataframe
    sample = df.sample(n=1)

    # Preprocess the sample in the same way as the training data
    sample_X = sample.drop(['Churn'], axis=1)
    sample_X_array = sample_X.values.reshape(1, -1)

    # Predict the income for the sample
    predicted_income = model.predict(sample_X_array)[0]

    # Store actual and predicted income in the results list
    results.append({
        "Sample": i+1,
        "Actual Income": sample['Churn'].values[0],
        "Predicted Income": predicted_income,
        "Difference": sample['Churn'].values[0] - predicted_income
    })
    
    # Create a DataFrame from the results
result_df = pd.DataFrame(results)
result_df








#linear reggression traing based on tenure and monthyl charges


df = pd.read_csv('WA_Fn-UseC_-Telco-Customer-Churn.csv')
df.info()

    # remove customerID 
df2 = df.drop(['customerID'], axis = 1)
df2.head(5)


# Categorical Data Encoding, Encode categorical variables, we use One-Hot Encoding for nominal variables and Label Encoding for ordinal variables.

#One-Hot Encoding: Gender, Partner, Dependents, MultipleLines, InternetService, OnlineSecurity, OnlineBackup, DeviceProtection, TechSupport, StreamingTV, StreamingMovies, PaperlessBilling, PaymentMethod, Churn

#Label Encoding: Contract
#fixing data types, e.g senior citinizen converted into yes and no

df2['TotalCharges'] =  pd.to_numeric(df2['TotalCharges'], errors='coerce')


# label encoding for binary columns
binary_cols = ['gender', 'SeniorCitizen', 'Dependents', 'Partner', 'Churn', 'PhoneService', 'PaperlessBilling']
le = LabelEncoder()
for col in binary_cols:
    df2[col] = le.fit_transform(df2[col])


# encoding remaining multiclass features using one hot encoding
multiclass_cols = ['OnlineSecurity', 'OnlineBackup', 'DeviceProtection', 'TechSupport', 'StreamingTV', 'StreamingMovies', 'MultipleLines', 'InternetService', 'Contract', 'PaymentMethod']

df3 = pd.get_dummies(df, columns=multiclass_cols, drop_first=True, dtype=np.int32)




# Label Encoding
from sklearn import preprocessing
label_encoder = preprocessing.LabelEncoder()

df.head()



print(df3.columns)
df4 = df[['Age',  'Survived',  'Fare', 'Pclass_2','Pclass_3', 'Embarked_Q', 'Embarked_S'] ]
df4.info()

# Define the features and the target variable
X = df4.drop(['Survived'], axis=1) #independent variables (predictors) 
y = df4['Survived'] #dependent variables (target


from collections import Counter
Counter(y)

## handling imbalanced data
from imblearn.combine import SMOTEENN
sm = SMOTEENN(random_state=0)
X, y = sm.fit_resample(X,y)


from collections import Counter
Counter(y)



# train test split
from sklearn.model_selection import train_test_split # split dataset

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Feature Scaling - Standardization 
 
standard_scaler = preprocessing.StandardScaler().fit(X_train)
X_train_standard = standard_scaler.transform(X_train)
X_test_standard = standard_scaler.transform(X_test)

     


#initialize the model
algorithm = LinearRegression() # here () different parameters can be use to optimize the model

model = algorithm.fit( X_train, y_train)



# predicting
y_pred = model.predict(X_test)
# Assuming y_pred contains probabilities
y_pred = np.where(y_pred > 0.5, 1, 0)


# evaluation
from sklearn.metrics import classification_report 
print(classification_report(y_test, y_pred))

#printing table for x_test and y_pred Assuming y_test and y_pred are arrays
data = {'y_test': y_test, 'y_pred': y_pred}
df_y_test_pred = pd.DataFrame(data)

# Display the table
print(df_y_test_pred)

# Save the table as a CSV file
df.to_csv('y_test_y_pred_table.csv', index=False)

# confusion_matrix, unction from scikit-learn is a tool used to evaluate the performance of a classification model by comparing the true labels (y_test) with the predicted labels (y_pred).
from sklearn.metrics import confusion_matrix
print(confusion_matrix(y_test, y_pred))


# Generate the confusion matrix
cm = confusion_matrix(y_test, y_pred)  # Ensure y_pred_binary contains binary predictions

# Extract counts from the confusion matrix
TN, FP, FN, TP = cm.ravel()  # Flatten the 2x2 matrix into individual values

# Data for the bar chart
categories = ['True Negatives (TN)', 'False Positives (FP)', 'False Negatives (FN)', 'True Positives (TP)']
values = [TN, FP, FN, TP]

# Plot the bar chart
plt.figure(figsize=(8, 6))
plt.bar(categories, values, color=['blue', 'orange', 'red', 'green'])
plt.title('Confusion Matrix Counts')
plt.ylabel('Count')
plt.xlabel('Category')
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.show()

#model Validation# Validate the model
y_val_pred = model.predict(X_test)
val_mse = mean_squared_error(y_test, y_val_pred) #mean squre error
accuracy = model.score(X_test, y_test) # accuracy (R-sqaured)

print(f"Validation Mean Squared Error: {val_mse}")

print(f"Validation Accuracy: {round(accuracy*100,2)}%")




#  Make predictions on the test set
y_test_pred = model.predict(X_test)

# Evaluate the model

test_mse = mean_squared_error(y_test, y_test_pred) #mean squre error
accuracy = model.score(X_test, y_test) # accuracy (R-sqaured)

print(f"Test Mean Squared Error: {test_mse}")

print(f"Testing Accuracy: {round(accuracy*100,2)}%")

# performance matrix
from sklearn.metrics import accuracy_score, f1_score, precision_score ,recall_score, roc_auc_score
accuracy = round(accuracy_score(y_test, y_pred),2)
f1_score = round(f1_score(y_test, y_pred),2)
precision = round(precision_score(y_test, y_pred),2)
recall = round(recall_score(y_test, y_pred),2)


#logis
from astropy.table import Table
dict1 = [{'accuracy': accuracy, 'f1_score': f1_score, 'precision': precision, 'recall': recall, }]
logis_matrix = Table(rows=dict1)
print(logis_matrix)

residuals = y_test - y_pred
plt.hist(residuals, bins=30, color='orange', alpha=0.7)
plt.title('Error Distribution')
plt.xlabel('Residuals')
plt.ylabel('Frequency')
plt.show()


#fpr new data
# NEW DATA: Predict on new data
# Example: new_data = [[value1, value2, value3, ...]] (ensure feature structure matches training data)
new_data_standard = standard_scaler.transform(new_data)  # Standardize new data
new_predictions = lm.predict(new_data_standard)  # Predict classes
new_probabilities = lm.predict_proba(new_data_standard)  # Predict probabilities

# Print predictions for new data
print("Predicted Classes:", new_predictions)
print("Predicted Probabilities:", new_probabilities)