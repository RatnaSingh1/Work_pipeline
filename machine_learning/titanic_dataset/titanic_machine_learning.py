# -*- coding: utf-8 -*-
"""
Created on Tue May 20 14:11:52 2025

@author: ratsa
"""

import pandas as pd
import numpy as np
from sklearn.preprocessing import OneHotEncoder, OrdinalEncoder, StandardScaler
from sklearn.compose import ColumnTransformer
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from imblearn.combine import SMOTEENN
from sklearn.linear_model import Perceptron
# Load dataset
df = pd.read_csv('train.csv')

# Remove unnecessary columns
df.drop(["Name", "Cabin", "SibSp", "Parch", "Ticket"], axis=1, inplace=True)

# Handle missing values
df['Age'].fillna(df['Age'].median(), inplace=True)
df['Embarked'].fillna('S', inplace=True)

# Convert age to categorical groups
def categorize_age(age):
    if age <= 4:
        return "Baby"
    elif age <= 12:
        return "Child"
    elif age <= 19:
        return "Teen"
    elif age <= 39:
        return "Adult"
    elif age <= 59:
        return "Middle Age Adult"
    else:
        return "Senior Adult"

df["Age_group"] = df["Age"].apply(categorize_age)
df.drop("Age", axis=1, inplace=True)

# Store PassengerId separately
passenger_ids = df["PassengerId"]
df.drop("PassengerId", axis=1, inplace=True)

# Encoding categorical features
Ohe = OneHotEncoder(drop='first', sparse_output=False)
OE = OrdinalEncoder(categories=[[1, 2, 3], ["Baby", "Child", "Teen", "Adult", "Middle Age Adult", "Senior Adult"]])
binary_encoder = OrdinalEncoder(categories=[["male", "female"]])

preprocessor = ColumnTransformer([
    ("onehot_encoding", Ohe, ["Embarked"]),
    ("ordinal_encoding", OE, ["Pclass", "Age_group"]),
    ("binary_encoding", binary_encoder, ["Sex"])
], remainder="passthrough")

df_encoded = preprocessor.fit_transform(df)
df_encoded = pd.DataFrame(df_encoded, columns=preprocessor.get_feature_names_out())

df_encoded.rename(columns={
    "remainder__Fare": "Fare", "binary_encoding__Sex": "Sex",
    "onehot_encoding__Embarked_Q": "Embarked_Q", "onehot_encoding__Embarked_S": "Embarked_S",
    "ordinal_encoding__Pclass": "Pclass", "ordinal_encoding__Age_group": "Age_group",
    "remainder__Survived": "Survived"
}, inplace=True)

# Feature scaling
scaler = StandardScaler()
df_encoded_scaled = scaler.fit_transform(df_encoded)
df_encoded = pd.DataFrame(df_encoded_scaled, columns=df_encoded.columns)

df_encoded["Survived"] = df_encoded["Survived"].astype(int)

# Splitting into features and target variable
X = df_encoded.drop(['Survived'], axis=1)
y = df_encoded['Survived']

# Handling imbalanced data
sm = SMOTEENN(random_state=0)
X_resampled, y_resampled = sm.fit_resample(X, y)

# Split dataset into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X_resampled, y_resampled, test_size=0.2, random_state=42)

# Train and evaluate models
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.neural_network import MLPClassifier
from sklearn.linear_model import Perceptron

models_titanic = {
    'Logistic Regression': LogisticRegression(max_iter=5000, random_state=42),
    'Random Forest': RandomForestClassifier(),
    'Gradient Boosting': GradientBoostingClassifier(),
    'Decision Tree': DecisionTreeClassifier(),
    'Support Vector Machine': SVC(),
    'K-Nearest Neighbors': KNeighborsClassifier(),
    'Naïve Bayes': GaussianNB(),
    'Perceptron': Perceptron(),
    'Multi-Layer Perceptron': MLPClassifier()
}

def evaluate_models(X_train, X_test, y_train, y_test):
    results = []
    predictions = {}
    
    for name, model in models_titanic.items():
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)
        acc = accuracy_score(y_test, y_pred)
        results.append((name, acc))
        predictions[name] = y_pred
    
    results.sort(key=lambda x: x[1], reverse=True)
    return results, predictions

results, predictions = evaluate_models(X_train, X_test, y_train, y_test)

# Best model selection
best_model_name = results[0][0]
best_model = models_titanic[best_model_name]
y_test_pred = best_model.predict(X_test)

# Restore PassengerId and add predictions
final_output = pd.DataFrame(X_test)
final_output['Survived'] = y_test_pred
final_output['PassengerId'] = passenger_ids.iloc[X_test.index]  # Ensure correct PassengerId mapping

# Save predictions with PassengerId
final_output.to_csv("titanic_predictions.csv", index=False)

# Print first few rows for manual inspection
print("\nActual vs. Predicted values for the best model:")
print(df_comparison.head())



predictions_df = pd.DataFrame(predictions)
print(predictions_df)


# Print model performance
print("Model Performance:")
for name, acc in results: print(f"{name}: {acc:.6f}")

# Create DataFrame to compare actual vs. predicted values for the best model
best_model_name = results[0][0]  # after soring (results.sort)Get the model with highest accuracy
y_test_pred = predictions[best_model_name]  # Get its predictions

df_comparison = pd.DataFrame({
    "Actual Value": y_test,
    "Predicted Value": y_test_pred
})

# Print first few rows for manual inspection
print("\nActual vs. Predicted values for the best model:")
print(df_comparison.head())

# Save to Excel for further inspection (optional)
df_comparison.to_excel("predictions.xlsx", index=False)
print("Saved predictions to predictions.xlsx")


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import OneHotEncoder, OrdinalEncoder, StandardScaler
from sklearn.compose import ColumnTransformer

# Load test dataset
new_data = pd.read_csv('test.csv')

# Store PassengerId separately before transformation
passenger_ids = new_data["PassengerId"]
new_data.drop("PassengerId", axis=1, inplace=True)

# Remove unnecessary columns
new_data.drop(["Name", "Cabin", "SibSp", "Parch", "Ticket"], axis=1, inplace=True)

# Handle missing values
new_data['Age'].fillna(new_data['Age'].median(), inplace=True)
new_data['Embarked'].fillna('S', inplace=True)
new_data['Fare'].fillna(new_data['Fare'].median(), inplace=True)

# Convert age to categorical groups
def categorize_age(age):
    if age <= 4:
        return "Baby"
    elif age <= 12:
        return "Child"
    elif age <= 19:
        return "Teen"
    elif age <= 39:
        return "Adult"
    elif age <= 59:
        return "Middle Age Adult"
    else:
        return "Senior Adult"

new_data["Age_group"] = new_data["Age"].apply(categorize_age)
new_data.drop("Age", axis=1, inplace=True)

# Encode categorical features
Ohe = OneHotEncoder(drop='first', sparse_output=False)
OE = OrdinalEncoder(categories=[[1, 2, 3], ["Baby", "Child", "Teen", "Adult", "Middle Age Adult", "Senior Adult"]])
binary_encoder = OrdinalEncoder(categories=[["male", "female"]])

preprocessor = ColumnTransformer([
    ("onehot_encoding", Ohe, ["Embarked"]),
    ("ordinal_encoding", OE, ["Pclass", "Age_group"]),
    ("binary_encoding", binary_encoder, ["Sex"])
], remainder="passthrough")

new_data_encoded = preprocessor.fit_transform(new_data)
new_data_encoded = pd.DataFrame(new_data_encoded, columns=preprocessor.get_feature_names_out())

# Rename columns for consistency
new_data_encoded.rename(columns={
    "remainder__Fare": "Fare", "binary_encoding__Sex": "Sex",
    "onehot_encoding__Embarked_Q": "Embarked_Q", "onehot_encoding__Embarked_S": "Embarked_S",
    "ordinal_encoding__Pclass": "Pclass", "ordinal_encoding__Age_group": "Age_group"
}, inplace=True)

# Apply scaling
scaler = StandardScaler()
new_data_scaled = scaler.fit_transform(new_data_encoded)
new_data_encoded = pd.DataFrame(new_data_scaled, columns=new_data_encoded.columns)

# Load trained model
from sklearn.ensemble import GradientBoostingClassifier
best_model = GradientBoostingClassifier()
best_model.fit(X_train, y_train)  # Ensure model is trained before predicting

# Make predictions
y_pred_new = best_model.predict(new_data_encoded)

# Restore PassengerId and add predictions
new_data_encoded['Predicted_Survived'] = y_pred_new
new_data_encoded['PassengerId'] = passenger_ids  # Reattach PassengerId

# Save final output with predictions
new_data_encoded.to_csv("titanic_predictions.csv", index=False)

# Plot survival distribution based on Fare
plt.figure(figsize=(8, 5))
sns.boxplot(data=new_data_encoded, x="Predicted_Survived", y="Fare")
plt.xlabel("Predicted Survived")
plt.ylabel("Fare")
plt.title("Fare Distribution by Survival")
plt.show()

# Categorize fares for better visualization
bins = [0, 100, 300, 600]
labels = ["Low (0-100)", "Medium (100-300)", "High (300-600)"]
new_data_encoded["Fare_Category"] = pd.cut(new_data_encoded["Fare"], bins=bins, labels=labels)

plt.figure(figsize=(8, 5))
sns.countplot(data=new_data_encoded, x="Fare_Category", hue="Predicted_Survived")
plt.xlabel("Fare Category")
plt.ylabel("Count")
plt.title("Survival Count by Fare Category")
plt.legend(title="Survived", labels=["Did Not Survive (0)", "Survived (1)"])
plt.show()
