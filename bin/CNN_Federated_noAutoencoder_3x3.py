#!/usr/bin/env python
import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px
import seaborn as sns
import tensorflow as tf
from tensorflow.keras import layers, models
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import mean_squared_error

def run_federated(args):
    studies = ["SDY569", "SDY797", "SDY1737"]
    standard_features = ["MIAA", "GAD65", "IA2IC", "ICA", "ZNT8", "8-12", "13-17", ">18", "Sex"]

    os.makedirs("models", exist_ok=True)
    os.makedirs("figures/pdf", exist_ok=True)
    os.makedirs("figures/html", exist_ok=True)

    # Customizable input train/test paths
    train_paths = {
        "SDY569": args.sdy569_train,
        "SDY797": args.sdy797_train,
        "SDY1737": args.sdy1737_train
    }
    test_paths = {
        "SDY569": args.sdy569_test,
        "SDY797": args.sdy797_test,
        "SDY1737": args.sdy1737_test
    }
    
    # 1. Load train/test data
    split_data = {}
    for study in studies:
        train_df = pd.read_csv(train_paths[study])
        test_df = pd.read_csv(test_paths[study])
        train_df = train_df.fillna(0.0)
        test_df = test_df.fillna(0.0)
        train_df["log_auc"] = np.log(train_df["C_Peptide_AUC_4Hrs"])
        test_df["log_auc"] = np.log(test_df["C_Peptide_AUC_4Hrs"])
        X_train = train_df[standard_features].values
        y_train = train_df["log_auc"].values
        X_test = test_df[standard_features].values
        y_test = test_df["log_auc"].values
        scaler = MinMaxScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        X_test_scaled = scaler.transform(X_test)
        X_train_cnn = X_train_scaled.reshape(-1, 3, 3, 1)
        X_test_cnn = X_test_scaled.reshape(-1, 3, 3, 1)
        split_data[study] = {
            "X_train": X_train_cnn,
            "y_train": y_train,
            "X_test": X_test_cnn,
            "y_test": y_test
        }

    # 2. Train local models & collect weights
    all_weights = []
    for study in studies:
        print(f"\nTraining local model for {study}")
        X_train = split_data[study]["X_train"]
        y_train = split_data[study]["y_train"]
        model = models.Sequential([
            layers.Input(shape=(3, 3, 1)),
            layers.Conv2D(16, kernel_size=(3, 3), activation='relu'),
            layers.Flatten(),
            layers.Dense(32, activation='relu'),
            layers.Dense(1)
        ])
        model.compile(optimizer='adam', loss='mse')
        model.fit(X_train, y_train, epochs=100, batch_size=8, verbose=0)
        all_weights.append(model.get_weights())

    # 3. Federated averaging (all studies)
    federated_weights = []
    for i in range(len(all_weights[0])):
        stacked = np.stack([w[i] for w in all_weights], axis=0)
        avg = np.mean(stacked, axis=0)
        federated_weights.append(avg)
    federated_model = models.Sequential([
        layers.Input(shape=(3, 3, 1)),
        layers.Conv2D(16, kernel_size=(3, 3), activation='relu'),
        layers.Flatten(),
        layers.Dense(32, activation='relu'),
        layers.Dense(1)
    ])
    federated_model.compile(optimizer='adam', loss='mse')
    federated_model.set_weights(federated_weights)
    federated_model.save("models/federated_CNN_3x3_noAutoencoder.keras")

    # 4. Evaluate federated model (all studies)
    for study in studies:
        print(f"\nEvaluating Federated Model on {study}")
        X_test = split_data[study]["X_test"]
        y_test = split_data[study]["y_test"]
        y_pred = federated_model.predict(X_test).flatten()
        mse = mean_squared_error(y_test, y_pred)
        median_mse = np.median((y_test - y_pred) ** 2)
        iqr = np.percentile((y_test - y_pred) ** 2, 75) - np.percentile((y_test - y_pred) ** 2, 25)
        print(f"  MSE:        {mse:.4f}")
        print(f"  Median MSE: {median_mse:.4f}")
        print(f"  IQR:        {iqr:.4f}")
        plt.figure(figsize=(6,6))
        plt.scatter(y_test, y_pred, alpha=0.8)
        plt.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'r--')
        plt.xlabel("True log(C-peptide AUC)")
        plt.ylabel("Predicted")
        plt.title(f"{study} — Federated CNN\nMSE={mse:.3f}, Med MSE={median_mse:.3f}, IQR={iqr:.3f}")
        plt.grid(True)
        plt.savefig(f"figures/pdf/{study}_federated_cnn.pdf")
        plt.show()
        df_plot = pd.DataFrame({"True": y_test, "Predicted": y_pred})
        fig = px.scatter(df_plot, x="True", y="Predicted", title=f"{study} — Federated CNN")
        fig.add_shape(
            type="line",
            x0=y_test.min(), y0=y_test.min(),
            x1=y_test.max(), y1=y_test.max(),
            line=dict(color="red", dash="dash")
        )
        fig.write_html(f"figures/html/{study}_federated_cnn.html")

    # 5. Federated averaging excluding SDY797
    studies_no_sdy797 = ["SDY569", "SDY1737"]
    all_weights_no_sdy797 = []
    for study in studies_no_sdy797:
        print(f"\nTraining local model for {study}")
        X_train = split_data[study]["X_train"]
        y_train = split_data[study]["y_train"]
        model = models.Sequential([
            layers.Input(shape=(3, 3, 1)),
            layers.Conv2D(16, kernel_size=(3, 3), activation='relu'),
            layers.Flatten(),
            layers.Dense(32, activation='relu'),
            layers.Dense(1)
        ])
        model.compile(optimizer='adam', loss='mse')
        model.fit(X_train, y_train, epochs=100, batch_size=8, verbose=0)
        all_weights_no_sdy797.append(model.get_weights())
    federated_weights_no_sdy797 = []
    for i in range(len(all_weights_no_sdy797[0])):
        stacked = np.stack([w[i] for w in all_weights_no_sdy797], axis=0)
        avg = np.mean(stacked, axis=0)
        federated_weights_no_sdy797.append(avg)
    federated_model_no_sdy797 = models.Sequential([
        layers.Input(shape=(3, 3, 1)),
        layers.Conv2D(16, kernel_size=(3, 3), activation='relu'),
        layers.Flatten(),
        layers.Dense(32, activation='relu'),
        layers.Dense(1)
    ])
    federated_model_no_sdy797.compile(optimizer='adam', loss='mse')
    federated_model_no_sdy797.set_weights(federated_weights_no_sdy797)
    federated_model_no_sdy797.save("models/federated_CNN_3x3_noAutoencoder_excluding_SDY797.keras")
    # Evaluate (excluding SDY797)
    for study in studies_no_sdy797:
        print(f"\nEvaluating Federated Model on {study} excluding SDY797")
        X_test = split_data[study]["X_test"]
        y_test = split_data[study]["y_test"]
        y_pred = federated_model_no_sdy797.predict(X_test).flatten()
        mse = mean_squared_error(y_test, y_pred)
        median_mse = np.median((y_test - y_pred) ** 2)
        iqr = np.percentile((y_test - y_pred) ** 2, 75) - np.percentile((y_test - y_pred) ** 2, 25)
        print(f"  MSE:        {mse:.4f}")
        print(f"  Median MSE: {median_mse:.4f}")
        print(f"  IQR:        {iqr:.4f}")
        plt.figure(figsize=(6,6))
        plt.scatter(y_test, y_pred, alpha=0.8)
        plt.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'r--')
        plt.xlabel("True log(C-peptide AUC)")
        plt.ylabel("Predicted")
        plt.title(f"{study} — Federated CNN excluding SDY797\nMSE={mse:.3f}, Med MSE={median_mse:.3f}, IQR={iqr:.3f}")
        plt.grid(True)
        plt.savefig(f"figures/pdf/{study}_federated_cnn_excluding_SDY797.pdf")
        plt.show()
        df_plot = pd.DataFrame({"True": y_test, "Predicted": y_pred})
        fig = px.scatter(df_plot, x="True", y="Predicted", title=f"{study} — Federated CNN excluding SDY797")
        fig.add_shape(
            type="line",
            x0=y_test.min(), y0=y_test.min(),
            x1=y_test.max(), y1=y_test.max(),
            line=dict(color="red", dash="dash")
        )
        fig.write_html(f"figures/html/{study}_federated_cnn_excluding_SDY797.html")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Federated CNN without Autoencoder 3x3 Filter")
    parser.add_argument("--working_dir", type=str, default=None, help="Working directory")
    parser.add_argument("--sdy569_train", type=str, default="data/cleaned/SDY569_train.csv")
    parser.add_argument("--sdy797_train", type=str, default="data/cleaned/SDY797_train.csv")
    parser.add_argument("--sdy1737_train", type=str, default="data/cleaned/SDY1737_train.csv")
    parser.add_argument("--sdy569_test", type=str, default="data/cleaned/SDY569_test.csv")
    parser.add_argument("--sdy797_test", type=str, default="data/cleaned/SDY797_test.csv")
    parser.add_argument("--sdy1737_test", type=str, default="data/cleaned/SDY1737_test.csv")
    args = parser.parse_args()
    if args.working_dir:
        os.chdir(args.working_dir)
    run_federated(args)