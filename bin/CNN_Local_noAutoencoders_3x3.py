#!/usr/bin/env python
import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler
from tensorflow.keras import layers, models

def main(args):
    # STEP 0: Set working directory if provided
    if args.working_dir:
        os.chdir(args.working_dir)
    
    # Feature Data Paths (raw tidy files)
    feature_paths = {
        "SDY569": args.sdy569_feature,
        "SDY797": args.sdy797_feature,
        "SDY1737": args.sdy1737_feature
    }
    # Ground Truth C-Peptide Labels
    cpeptide_paths = {
        "SDY569": args.sdy569_cpeptide,
        "SDY797": args.sdy797_cpeptide,
        "SDY1737": args.sdy1737_cpeptide
    }
    # CNN Output Features (Used for Ridge Regression)
    cnn_test_feature_paths = {
        "SDY569": args.sdy569_cnn_test,
        "SDY797": args.sdy797_cnn_test,
        "SDY1737": args.sdy1737_cnn_test
    }
    
    # Standard features
    standard_features = [
        "MIAA", "GAD65", "IA2IC", "ICA", "ZNT8",
        "8-12", "13-17", ">18", "Sex"
    ]
    cleaned_data = {}
    for study in feature_paths:
        print(f"\n=== Preprocessing {study} ===")
        df_feat = pd.read_csv(feature_paths[study])
        df_cpep = pd.read_csv(cpeptide_paths[study])
        if "Accession" in df_feat.columns:
            df_feat = df_feat.rename(columns={"Accession": "Subject_ID"})
        # One-hot encode Age_Group
        df_feat = pd.get_dummies(df_feat, columns=["Age_Group"])
        # Pivot
        df_wide = df_feat.pivot_table(
            index=["Subject_ID", "Sex"],
            columns="Property",
            values="Value"
        ).reset_index()
        # Map sex to binary
        df_wide["Sex"] = df_wide["Sex"].map({"Male": 0, "Female": 1})
        # Manually add expected age group indicators
        for col in ["Age_Group_8-12", "Age_Group_13-17", "Age_Group_>18"]:
            if col not in df_feat.columns:
                df_feat[col] = 0
        df_age = df_feat.groupby("Subject_ID")[["Age_Group_8-12", "Age_Group_13-17", "Age_Group_>18"]].max().reset_index()
        df_wide = pd.merge(df_wide, df_age, on="Subject_ID", how="left")
        df_wide = df_wide.rename(columns={
            "Age_Group_8-12": "8-12",
            "Age_Group_13-17": "13-17",
            "Age_Group_>18": ">18"
        })
        df = pd.merge(df_wide, df_cpep, on="Subject_ID", how="inner")
        for feat in standard_features:
            if feat not in df.columns:
                df[feat] = 0.0
        df = df[["Subject_ID"] + standard_features + ["C_Peptide_AUC_4Hrs"]]
        print(df.head(3))
        cleaned_data[study] = df
    
    # Ensure output directories exist
    os.makedirs("figures/pdf", exist_ok=True)
    os.makedirs("figures/html", exist_ok=True)
    os.makedirs("models", exist_ok=True)
    os.makedirs("data/cleaned", exist_ok=True)
    
    # Local Training Loop
    for study, df in cleaned_data.items():
        print(f"\n=== Training Local CNN for {study} ===")
        df = df.fillna(0.0)
        df["log_auc"] = np.log(df["C_Peptide_AUC_4Hrs"])
        X = df[standard_features].values
        y = df["log_auc"].values

        print(X)
        print(y)

        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.5, random_state=42)
        # Save cleaned full dataset
        df.to_csv(f"data/cleaned/{study}_cleaned.csv", index=False)
        subject_ids = df["Subject_ID"].values
        train_ids = subject_ids[:len(X_train)]
        test_ids = subject_ids[len(X_train):]
        df_train = df[df["Subject_ID"].isin(train_ids)]
        df_test = df[df["Subject_ID"].isin(test_ids)]
        df_train.to_csv(f"data/cleaned/{study}_train.csv", index=False)
        df_test.to_csv(f"data/cleaned/{study}_test.csv", index=False)

        scaler = MinMaxScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled  = scaler.transform(X_test)
        X_train_cnn = X_train_scaled.reshape(-1, 3, 3, 1)
        X_test_cnn  = X_test_scaled.reshape(-1, 3, 3, 1)

        # Define CNN
        model = models.Sequential([
            layers.Input(shape=(3, 3, 1)),
            layers.Conv2D(16, kernel_size=(3, 3), activation='relu'),
            layers.Flatten(),
            layers.Dense(32, activation='relu'),
            layers.Dense(1)
        ])
        model.compile(optimizer='adam', loss='mse')
        model.fit(X_train_cnn, y_train, epochs=100, batch_size=8, verbose=0)
        y_pred = model.predict(X_test_cnn).flatten()
        model.save(f"models/{study}_local_cnn.keras")

        squared_errors = (y_test - y_pred) ** 2
        avg_mse = np.mean(squared_errors)
        median_mse = np.median(squared_errors)
        iqr = np.percentile(squared_errors, 75) - np.percentile(squared_errors, 25)
        print(f"  Avg MSE:    {avg_mse:.4f}")
        print(f"  Median MSE: {median_mse:.4f}")
        print(f"  IQR:        {iqr:.4f}")

        plt.figure(figsize=(6, 6))
        plt.scatter(y_test, y_pred, alpha=0.8)
        plt.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'r--')
        plt.xlabel("True log(C-peptide AUC)")
        plt.ylabel("Predicted")
        plt.title(f"{study} — Local CNN\nAvg MSE={avg_mse:.3f}, Med MSE={median_mse:.3f}, IQR={iqr:.3f}")
        plt.grid(True)
        plt.savefig(f"figures/pdf/{study}_local_cnn.pdf")
        plt.show()

        # Interactive plot with plotly
        df_plot = pd.DataFrame({
            "True log(AUC)": y_test,
            "Predicted": y_pred
        })
        fig = px.scatter(
            df_plot, x="True log(AUC)", y="Predicted",
            title=f"{study} — Local CNN (Interactive)",
            labels={"x": "True log(C‑peptide AUC)", "y": "Predicted"},
            width=600, height=600
        )
        fig.add_shape(
            type="line",
            x0=y_test.min(), y0=y_test.min(),
            x1=y_test.max(), y1=y_test.max(),
            line=dict(color="red", dash="dash")
        )
        fig.write_html(f"figures/html/{study}_local_cnn.html")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="CNN for C Peptide Prediction (No Autoencoders)")
    parser.add_argument("--working_dir", type=str, default=None, help="Working directory to set at start of script")
    parser.add_argument("--sdy569_feature", type=str, default="data/SDY569_tidy.csv", help="Feature file for SDY569")
    parser.add_argument("--sdy797_feature", type=str, default="data/SDY797_tidy.csv", help="Feature file for SDY797")
    parser.add_argument("--sdy1737_feature", type=str, default="data/SDY1737_tidy.csv", help="Feature file for SDY1737")
    parser.add_argument("--sdy569_cpeptide", type=str, default="data/SDY569_cpeptide_auc_tidy.csv", help="C-peptide label file for SDY569")
    parser.add_argument("--sdy797_cpeptide", type=str, default="data/SDY797_cpeptide_auc_tidy.csv", help="C-peptide label file for SDY797")
    parser.add_argument("--sdy1737_cpeptide", type=str, default="data/SDY1737_cpeptide_auc_tidy.csv", help="C-peptide label file for SDY1737")
    parser.add_argument("--sdy569_cnn_test", type=str, default="test_data/SDY569_3x3_test.csv", help="CNN test feature file for SDY569")
    parser.add_argument("--sdy797_cnn_test", type=str, default="test_data/SDY797_3x3_test.csv", help="CNN test feature file for SDY797")
    parser.add_argument("--sdy1737_cnn_test", type=str, default="test_data/SDY1737_3x3_test.csv", help="CNN test feature file for SDY1737")
    args = parser.parse_args()
    main(args)