"""
Train a new LightGBM model for predicting bike sharing demand
"""

import pandas as pd
from lightgbm import LGBMRegressor
import mlflow
import os
from typing import Dict

def pull_data(dataset_url: str) -> tuple:
    """
    Download the data set from a given url
    Args: 
        dataset_url: dataset url
    Returns:
        A Pandas DataFrame of the dataset
    """
    input_df = pd.read_csv(dataset_url)
    return input_df

def preprocess(input_df: pd.DataFrame) -> (pd.DataFrame, pd.DataFrame):
    """
    Preprocess the data and split them into a training and a test dataset
    Args:
        input_df: The DataFrame of the whole dataset
    Returns:
        Two DataFrames, one for training and another for testing dataset
    """
    input_df["datetime"] = pd.to_datetime(input_df["datetime"])

    #create hour, day and month variables from datetime column
    input_df["hour"] = input_df["datetime"].dt.hour
    input_df["day"] = input_df["datetime"].dt.day
    input_df["month"] = input_df["datetime"].dt.month

    #drop datetime, casual and registered columns
    input_df.drop(["datetime", "casual", "registered"], axis=1, inplace=True)
    
    #split the original dataset into a training and a test dataset
    horizon = 24 * 7
    train, test = input_df.iloc[:-horizon,:], input_df.iloc[-horizon:,:]
    return train, test

# mlflow configuration
MLFLOW_S3_ENDPOINT_URL = "http://mlflow-minio.local"
MLFLOW_TRACKING_URI = "http://mlflow-server.local"
AWS_ACCESS_KEY_ID = "minioadmin"
AWS_SECRET_ACCESS_KEY = "minioadmin"
mlflow_experiment_name = "week1-lgbm-bike"

os.environ["MLFLOW_S3_ENDPOINT_URL"] = MLFLOW_S3_ENDPOINT_URL
os.environ["AWS_ACCESS_KEY_ID"] = AWS_ACCESS_KEY_ID
os.environ["AWS_SECRET_ACCESS_KEY"] = AWS_SECRET_ACCESS_KEY
mlflow.set_tracking_uri(MLFLOW_TRACKING_URI)
mlflow.set_experiment(mlflow_experiment_name)

# dataset URL
dataset_url = "https://raw.githubusercontent.com/yumoL/mlops_eng_course_datasets/master/intro/bike-demanding/train_full.csv"

# name of the target column
target = "count"

def train(params: Dict):
    """
    Train a LightGBM regression model, record the training and upload the model to MLflow service
    Args:
        params: Parameters passed to the LightGBM model
    """
    input_df = pull_data(dataset_url)
    train, _ = preprocess(input_df)

    # Prepare train_x and train_y
    train_x = train.drop([target], axis=1)
    train_y = train[[target]]

    with mlflow.start_run() as run:
        model = LGBMRegressor(**params)
        model.fit(train_x, train_y)

        for hyperparam_name, value in params.items():
            mlflow.log_param(hyperparam_name, value)

        model_name = "lgbm-bike"
        mlflow.lightgbm.log_model(
            lgb_model=model,
            artifact_path=model_name,
            registered_model_name="Week1LgbmBikeDemand"
        )

        print(f"The trained model is located at {mlflow.get_artifact_uri(artifact_path=model_name)}")


