#!/usr/bin/env python
"""
Generic single-study local CNN training script for endotype analysis.
This script processes one study at a time and prepares outputs for federated aggregation.
"""
import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler
from tensorflow.keras import layers, models

def main():
    parser = argparse.ArgumentParser(description='Local CNN training for single study')
    parser.add_argument('--study_id', required=True, help='Study identifier (e.g., SDY569)')
    parser.add_argument('--study_feature', required=True, help='Path to tidy feature CSV file')
    parser.add_argument('--study_cpeptide', required=True, help='Path to C-peptide AUC CSV file')
    parser.add_argument('--output_dir', default='./', help='Output directory')
    parser.add_argument('--working_dir', help='Working directory to change to')
    
    args = parser.parse_args()
    
    # STEP 0: Set working directory if provided
    if args.working_dir:
        os.chdir(args.working_dir)
    
    print(f"\n=== Processing {args.study_id} ===")
    
    # Create output directories
    os.makedirs("models", exist_ok=True)
    os.makedirs("figures/pdf", exist_ok=True)
    os.makedirs("figures/html", exist_ok=True)
    os.makedirs("test_data", exist_ok=True)
    
    # Standard features
    standard_features = [
        "MIAA", "GAD65", "IA2IC", "ICA", "ZNT8",
        "8-12", "13-17", ">18", "Sex"
    ]
    
    # Load data
    df_feat = pd.read_csv(args.study_feature)
    df_cpep = pd.read_csv(args.study_cpeptide)
    
    # Standardize column names
    if "Accession" in df_feat.columns:
        df_feat = df_feat.rename(columns={"Accession": "Subject_ID"})
    
    # One-hot encode Age_Group
    df_feat = pd.get_dummies(df_feat, columns=["Age_Group"])
    
    # Pivot the data from long to wide format
    df_wide = df_feat.pivot_table(
        index=["Subject_ID", "Sex"],
        columns="Property",
        values="Value"
    ).reset_index()
    
    # Map sex to binary
    df_wide["Sex"] = df_wide["Sex"].map({"Male": 0, "Female": 1})
    
    # Manually add expected age group indicators if not present
    for col in ["Age_Group_8-12", "Age_Group_13-17", "Age_Group_>18"]:
        if col not in df_feat.columns:
            df_feat[col] = 0
    
    # Get age group information
    df_age = df_feat.groupby("Subject_ID")[["Age_Group_8-12", "Age_Group_13-17", "Age_Group_>18"]].max().reset_index()
    df_wide = pd.merge(df_wide, df_age, on="Subject_ID", how="left")
    
    # Rename age group columns to match expected format
    df_wide = df_wide.rename(columns={
        "Age_Group_8-12": "8-12",
        "Age_Group_13-17": "13-17",
        "Age_Group_>18": ">18"
    })
    
    # Merge with C-peptide data
    df_merged = pd.merge(df_wide, df_cpep, on="Subject_ID", how="inner")
    
    # Ensure all expected features are present
    for feat in standard_features:
        if feat not in df_merged.columns:
            df_merged[feat] = 0.0
    
    # Select only the required columns
    df_merged = df_merged[["Subject_ID"] + standard_features + ["C_Peptide_AUC_4Hrs"]]
    df_merged = df_merged.fillna(0.0)
    
    print(f"Data shape after preprocessing: {df_merged.shape}")
    print(f"Available columns: {list(df_merged.columns)}")
    print(df_merged.head(3))
    
    # Calculate log AUC
    df_merged["log_auc"] = np.log(df_merged["C_Peptide_AUC_4Hrs"])
    
    # Prepare features and target
    X = df_merged[standard_features].values
    y = df_merged["log_auc"].values
    
    # Scale features
    scaler = MinMaxScaler()
    X_scaled = scaler.fit_transform(X)
    
    # Train/test split
    X_train, X_test, y_train, y_test = train_test_split(
        X_scaled, y, test_size=0.2, random_state=42
    )
    
    # Reshape for CNN (3x3 grid)
    X_train_cnn = X_train.reshape(-1, 3, 3, 1)
    X_test_cnn = X_test.reshape(-1, 3, 3, 1)
    
    # Build CNN model
    model = models.Sequential([
        layers.Conv2D(32, (2, 2), activation='relu', input_shape=(3, 3, 1)),
        layers.MaxPooling2D((2, 2)),
        layers.Flatten(),
        layers.Dense(64, activation='relu'),
        layers.Dropout(0.2),
        layers.Dense(32, activation='relu'),
        layers.Dense(1)
    ])
    
    model.compile(optimizer='adam', loss='mse', metrics=['mae'])
    
    # Train the model
    history = model.fit(
        X_train_cnn, y_train,
        epochs=100,
        batch_size=32,
        validation_data=(X_test_cnn, y_test),
        verbose=1
    )
    
    # Save model
    model.save(f"models/{args.study_id}_local_cnn_model.h5")
    
    # Predictions
    train_pred = model.predict(X_train_cnn)
    test_pred = model.predict(X_test_cnn)
    
    # Save training and test data with predictions for federated aggregation
    train_df = pd.DataFrame({
        'Subject_ID': range(len(X_train)),
        'study_id': args.study_id,
        'actual': y_train,
        'predicted': train_pred.flatten(),
        'dataset': 'train'
    })
    
    test_df = pd.DataFrame({
        'Subject_ID': range(len(X_test)),
        'study_id': args.study_id,
        'actual': y_test,
        'predicted': test_pred.flatten(),
        'dataset': 'test'
    })
    
    # Combine and save results
    results_df = pd.concat([train_df, test_df], ignore_index=True)
    results_df.to_csv(f"{args.study_id}_local_results.csv", index=False)
    
    # Save cleaned train/test splits for potential federated use
    train_features_df = pd.DataFrame(X_train, columns=[f'feature_{i}' for i in range(X_train.shape[1])])
    train_features_df['C_Peptide_AUC_4Hrs'] = np.exp(y_train)
    train_features_df.to_csv(f"{args.study_id}_train.csv", index=False)
    
    test_features_df = pd.DataFrame(X_test, columns=[f'feature_{i}' for i in range(X_test.shape[1])])
    test_features_df['C_Peptide_AUC_4Hrs'] = np.exp(y_test)
    test_features_df.to_csv(f"{args.study_id}_test.csv", index=False)
    
    # Create visualization
    plt.figure(figsize=(10, 6))
    plt.plot(history.history['loss'], label='Training Loss')
    plt.plot(history.history['val_loss'], label='Validation Loss')
    plt.xlabel('Epochs')
    plt.ylabel('Loss')
    plt.title(f'{args.study_id} - Local CNN Training Loss')
    plt.legend()
    plt.savefig(f"figures/pdf/{args.study_id}_training_loss.pdf")
    plt.savefig(f"figures/{args.study_id}_training_loss.png")
    plt.close()
    
    # Scatter plot of predictions vs actual
    plt.figure(figsize=(8, 8))
    plt.scatter(y_test, test_pred, alpha=0.6)
    plt.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'r--', lw=2)
    plt.xlabel('Actual')
    plt.ylabel('Predicted')
    plt.title(f'{args.study_id} - Test Set Predictions')
    plt.savefig(f"figures/pdf/{args.study_id}_predictions.pdf")
    plt.savefig(f"figures/{args.study_id}_predictions.png")
    plt.close()
    
    # Calculate and print metrics
    from sklearn.metrics import mean_squared_error, r2_score
    train_mse = mean_squared_error(y_train, train_pred)
    test_mse = mean_squared_error(y_test, test_pred)
    test_r2 = r2_score(y_test, test_pred)
    
    print(f"\n=== {args.study_id} Results ===")
    print(f"Training MSE: {train_mse:.4f}")
    print(f"Test MSE: {test_mse:.4f}")
    print(f"Test R²: {test_r2:.4f}")
    print(f"Number of training samples: {len(X_train)}")
    print(f"Number of test samples: {len(X_test)}")
    
    # Save metrics
    metrics_df = pd.DataFrame({
        'study_id': [args.study_id],
        'train_mse': [train_mse],
        'test_mse': [test_mse],
        'test_r2': [test_r2],
        'n_train': [len(X_train)],
        'n_test': [len(X_test)]
    })
    metrics_df.to_csv(f"{args.study_id}_metrics.csv", index=False)

if __name__ == "__main__":
    main()