import pickle
from pathlib import Path
from typing import Tuple, List, Optional
import logging
import pandas as pd

import great_expectations as gx
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer
from sklearn.preprocessing import TargetEncoder

from thefuzz import fuzz, process

from great_expectations.data_context import FileDataContext
from great_expectations.datasource.fluent import BatchRequest
from great_expectations.checkpoint.types.checkpoint_result import CheckpointResult


# Random seed for making the assignments reproducible
RANDOM_SEED = 42

### START CODE HERE
# Copy the following functions here: file_reader, dataframe_merger, drop_futile_columns, correct_distance_unit, 
# string_transformer, typo_fixer, data_extraction, batch_creator, create_checkpoint, separate_X_and_y, impute_missing, 
# datetime_decomposer, condition_encoder, target_encode, store_features, store_targets, store_encoder, drop_unimportant_features, 
# feature_engineering_pipeline.
### END CODE HERE

def etl(path: Path, 
        gx_context_root_dir: Path, 
        gx_datasource_name: str, 
        gx_checkpoint_name: str,
        gx_expectation_suite_name: str,
        gx_run_name: str,
        feature_store_path: Path, 
        feature_file_name: str,
        encoder_file_name: str, 
        target_file_name: Optional[str]=None, 
        fit_encoder: bool=False, 
        targets_included: bool=True
    ):
    """
    This function loads, merges, cleans, and validates the specified data, extract features, and save the features in the feature store.
    Args:
        path (Path): Path to the folder where the files "deals.csv" and "housing_info.json" exist
        gx_context_root_dir (Path): The directory containing the Great Expectations configs.
        gx_datasource_name (str): Great Expectations data source name
        gx_checkpoint_name (str): Name of the Great Expectations checkpoint that runs the validation
        gx_expectation_suite_name (str): Name of the Expectation Suite used to validate the data
        gx_run_name (str): Name of the validation running
        feature_store_path (Path): Path of the feature store
        feature_file_name (str): Filename for the stored features
        encoder_file_name (str): Filename for the stored target encoder
        target_file_name (None|str): Filename for the stored targets
        fit_encoder (bool): Whether a new target encoder should be fitted. If False, uses a previously stored encoder
        targets_included (bool): If True, df has all of the housing data including targets. If False, df has only the features.
    """
    correct_condition_values=['poor', 'tolerable', 'satisfactory', 'good', 'excellent']
    ### START CODE HERE
    
    ### END CODE HERE


