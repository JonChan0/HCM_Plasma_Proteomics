import numpy as np
import pandas as pd
from sklearn.impute import KNNImputer
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.feature_selection import SelectPercentile, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier
from sklearn.svm import SVC
from sklearn.compose import ColumnTransformer
from sklearn.model_selection import GridSearchCV, train_test_split, StratifiedKFold
import joblib
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import accuracy_score, roc_curve, auc, precision_recall_curve, confusion_matrix, f1_score
import os
import wandb
from sklearn.pipeline import Pipeline
import time
import sklearn
from sklearn.cross_decomposition import PLSRegression
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.base import BaseEstimator, TransformerMixin

sklearn.set_config(transform_output="pandas")

class PLSComponentTransformer(BaseEstimator, TransformerMixin):
        def __init__(self, n_components=2):
            self.n_components = n_components
            self.pls = PLSRegression(n_components=self.n_components)
        
        def fit(self, X, y=None):
            self.pls.fit(X, y)
            return self
        
        def transform(self, X):
            return self.pls.transform(X)

# Function to define the model and parameter grid based on user input
def get_model_and_params(model_type):
    if model_type in ['logistic_regression', 'logistic_regression_no_fs']:
        model = LogisticRegression(solver = 'liblinear')
        param_grid = {
            'classifier__C': [0.1, 1, 10],
            'classifier__class_weight': ['balanced']
        }
    elif model_type in ['random_forest', 'random_forest_no_fs']:
        model = RandomForestClassifier(random_state=42)
        param_grid = {
            'classifier__n_estimators': [100, 200],
            'classifier__max_depth': [None, 10, 20],
            'classifier__class_weight': ['balanced']
        }
    elif model_type in ['xgboost', 'xgboost_no_fs']:
        model = XGBClassifier(use_label_encoder=False, eval_metric='logloss', n_jobs=1, random_state=42)
        param_grid = {
            'classifier__n_estimators': [100, 200],
            'classifier__learning_rate': [0.01, 0.1, 0.2],
            'classifier__scale_pos_weight': [1, 10, 100, 1000]
        }
    elif model_type in ['svm', 'svm_no_fs']:
        model = SVC(probability=True,random_state=42)
        param_grid = {
            'classifier__C': [0.1, 1, 10],
            'classifier__kernel': ['linear', 'rbf', 'poly'],
            'classifier__class_weight': ['balanced']
        }
    elif model_type == 'l1_logistic_regression':
        model = LogisticRegression(penalty='l1', solver='liblinear',random_state=42)
        param_grid = {
            'classifier__C': [0.1, 1, 10],
            'classifier__class_weight': ['balanced']
        }
    elif model_type == 'elastic_net_logistic_regression':
        model = LogisticRegression(penalty='elasticnet', solver='saga', l1_ratio=0.5, random_state=42)
        param_grid = {
            'classifier__C': [0.1, 1, 10],
            'classifier__class_weight': ['balanced'],
            'classifier__l1_ratio': [0.1, 0.5, 0.9]
        }
    elif model_type == 'spls_lda':
        model = Pipeline([
            ('pls', PLSComponentTransformer()),  # Use the custom transformer
            ('lda', LinearDiscriminantAnalysis())
        ])
        param_grid = {
            'classifier__pls__n_components': [2, 5, 10, 30, 50, 100, 200]
        }
    else:
        raise ValueError("Unsupported model type. Choose from 'logistic_regression', 'random_forest', 'xgboost', 'svm', 'l1_logistic_regression', 'spls_lda', 'elastic_net_logistic_regression'.")
    
    return model, param_grid

# Generalized function to create and train the model pipeline with class weights
def train_model(X_train, y_train, model, param_grid, model_name, model_output_folder, X_train_data_path, feature_selection='False', features_to_bypass_fs=[], features_to_select_fs=[]):
    start_time = time.time()

    # Determine feature names for quantitative (numeric) features only in the original X_train prior to filtering for those in features_to_bypass_fs and features_to_select_fs
    quantitative_feature_names = X_train.select_dtypes(include=['int64', 'float64']).columns.tolist()

    # Determine feature names for categorical features only
    categorical_feature_names = X_train.select_dtypes(include=['object']).columns.tolist()

    # Determine feature names for boolean features
    boolean_feature_names = X_train.select_dtypes(include=['bool']).columns.tolist()

    print(quantitative_feature_names)
    print(categorical_feature_names)
    print(boolean_feature_names)

    features_to_select_names = [col for col in quantitative_feature_names if col in features_to_select_fs]
    features_to_bypass_names = [col for col in quantitative_feature_names if col in features_to_bypass_fs]

    ##############################################################################
    #From the X_train, only select the features that are in the features_to_bypass_fs and features_to_select_fs lists as well as in the categorical and boolean features
    X_train = X_train[features_to_bypass_names +  categorical_feature_names + boolean_feature_names + features_to_select_names ]
    print(X_train.columns)

    quantitative_feature_names = X_train.select_dtypes(include=['int64', 'float64']).columns.tolist() #Redefine the quantitative feature names to ensure that the correct ones are passed to KNN impute
    ###########################################################################################################
      # Check if model type is random forest or xgboost
    if model_name in ['random_forest', 'xgboost'] and feature_selection == 'True':
        bypass_pipeline = Pipeline(steps=[('pass','passthrough')])  # No scaling
        selection_pipeline = Pipeline(steps=[
            ('feature_selection', SelectPercentile(score_func=f_classif, percentile=1))  # Feature selection only
        ])
    elif model_name in ['random_forest', 'xgboost'] and feature_selection == 'False':
        bypass_pipeline = Pipeline(steps=[('pass','passthrough')]) #No scaling
        selection_pipeline = Pipeline(steps=[('pass','passthrough')]) #No feature selection
    
    else: #For other model types beyond tree-based
        # Preprocessing for bypassed quantitative features (scaling)
        bypass_pipeline = Pipeline(steps=[
            ('scaler', StandardScaler())  
        ])

        if feature_selection == 'False' or model_name in ['l1_logistic_regression', 'elastic_net_logistic_regression']:
            # Preprocessing and feature selection for selected quantitative features
            selection_pipeline = Pipeline(steps=[
                ('scaler', StandardScaler()) 
            ])
        else:
              # Preprocessing and feature selection for selected quantitative features
            selection_pipeline = Pipeline(steps=[
                ('scaler', StandardScaler()),
                ('feature_selection', SelectPercentile(score_func=f_classif, percentile=1))   #Feature selection = SelectPercentile
            ])
    ###########################################################################################################
    # Apply KNNImputer to all quantitative features and OneHotEncoder to all categorical and boolean features
    preprocessor_before_split = ColumnTransformer(
        transformers=[
            ('imputer', KNNImputer(n_neighbors=5), quantitative_feature_names),
            ('encoder', OneHotEncoder(handle_unknown='ignore', sparse_output=False), categorical_feature_names),
            ('boolean_pass', Pipeline(steps=[('pass','passthrough')]),boolean_feature_names)
        ],verbose_feature_names_out=False
    )

    # Combine the pipelines for bypass and for select features
    quantitative_pipeline = ColumnTransformer(
        transformers=[
            ('bypass', bypass_pipeline, features_to_bypass_names),
            ('select', selection_pipeline, features_to_select_names),
            ('boolean_pass', Pipeline(steps=[('pass','passthrough')]),boolean_feature_names)
        ],verbose_feature_names_out=False, remainder='passthrough' #This remainder = one-hot encoded categorical variables
    )
    ###########################################################################################################
    # Define the final pipeline with KNNImputer applied before splitting
    pipeline = Pipeline(steps=[
        ('imputer_preprocessor', preprocessor_before_split),  # Impute all quantitative features and scale categorical features
        ('feature_preprocessor', quantitative_pipeline),  # Process bypass and selected features + pass through categorical features
        ('classifier', model)  # Replace `model` with your classifier
    ])

    X_train_preprocessed = pipeline.named_steps['imputer_preprocessor'].fit_transform(X_train)
    X_train_preprocessed = pipeline.named_steps['feature_preprocessor'].fit_transform(X_train_preprocessed, y_train)
    base_folder = os.path.dirname(X_train_data_path)
    preprocessed_data_path = os.path.join(base_folder, f'X_train_preprocessed_{model_name}.csv')

    pd.DataFrame(X_train_preprocessed, columns=pipeline.named_steps['feature_preprocessor'].get_feature_names_out()).to_csv(preprocessed_data_path, index=False)
    print(f"Preprocessed X_train data saved to {preprocessed_data_path}")
    wandb.log({"status": "Saved preprocessed training data"})

    print("Starting model training...")
    wandb.log({"status": "Starting model training"})

    # Perform 5-fold cross-validation with parallelization and verbose output
    grid_search_start_time = time.time()

    grid_search = GridSearchCV(pipeline, param_grid, cv=StratifiedKFold(n_splits=5, shuffle=True, random_state=42), 
                               scoring=['average_precision','roc_auc', 'f1', 'precision','recall','balanced_accuracy','neg_brier_score'],
                               refit='roc_auc',
                               n_jobs=-1, verbose=10)
    grid_search.fit(X_train, y_train)
    grid_search_end_time = time.time()

    print("Model training completed.")
    wandb.log({"status": "Model training completed"})

    # Log the best parameters and best score to wandb
    wandb.log({"best_params": grid_search.best_params_, "best_score": grid_search.best_score_})

    # Print the best parameters and best score
    print(f"Best parameters found for {model.__class__.__name__}: ", grid_search.best_params_)
    print(f"Best mean cross-validation AUC for {model.__class__.__name__}: ", grid_search.best_score_)

    #Save the cv_results of the grid search object as well
    cv_results_path = os.path.join(model_output_folder, f'{model_name}_cv_results.csv')
    pd.DataFrame(grid_search.cv_results_).to_csv(cv_results_path, index=False)
    print(f"CV results saved to {cv_results_path}")

    # Save the best estimator
    model_path = os.path.join(model_output_folder, f'{model_name}_best_model.pkl')
    joblib.dump(grid_search.best_estimator_, model_path)

    # Log the model to wandb
    # wandb.save(f'{model_name}_best_model.pkl')

    # Log timing information
    total_time = time.time() - start_time
    grid_search_time = grid_search_end_time - grid_search_start_time
    wandb.log({"total_time": total_time, "grid_search_time": grid_search_time})
    print(f"Total time: {total_time:.2f} seconds")
    print(f"Grid search time: {grid_search_time:.2f} seconds")

    # Log summary metrics
    wandb.summary['best_score'] = grid_search.best_score_
    wandb.summary['total_time'] = total_time
    wandb.summary['grid_search_time'] = grid_search_time

    # Return the best estimator
    return grid_search.best_estimator_

# Function to plot metrics
def plot_metrics(model, X, y, dataset_name, model_name, plot_output_folder):

    #If plot_output_folder does not exist, create it:
    if not os.path.exists(plot_output_folder):
        os.makedirs(plot_output_folder)
    
    y_pred = model.predict(X)
    y_pred_proba = model.predict_proba(X)[:, 1]

    # Calculate accuracy
    accuracy = accuracy_score(y, y_pred)
    wandb.log({f"{dataset_name}_accuracy": accuracy})
    print(f"{dataset_name} set accuracy: ", accuracy)

    # Plot ROC curve
    fpr, tpr, _ = roc_curve(y, y_pred_proba)
    roc_auc = auc(fpr, tpr)
    plt.figure()
    plt.plot(fpr, tpr, color='darkorange', lw=2, label='ROC curve (area = %0.2f)' % roc_auc)
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(f'Receiver Operating Characteristic - {dataset_name}')
    plt.legend(loc="lower right")
    roc_curve_path = os.path.join(plot_output_folder, f'{model_name}_roc_curve_{dataset_name}.png')
    plt.savefig(roc_curve_path)
    plt.close()
    wandb.log({f"{dataset_name}_roc_curve": wandb.Image(roc_curve_path)})

    # Plot Precision-Recall curve
    precision, recall, _ = precision_recall_curve(y, y_pred_proba)
    plt.figure()
    plt.plot(recall, precision, color='blue', lw=2)
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title(f'Precision-Recall curve - {dataset_name}')
    pr_curve_path = os.path.join(plot_output_folder, f'{model_name}_precision_recall_curve_{dataset_name}.png')
    plt.savefig(pr_curve_path)
    plt.close()
    wandb.log({f"{dataset_name}_precision_recall_curve": wandb.Image(pr_curve_path)})

    # Confusion matrix and F1 score at different thresholds
    thresholds = [0.3, 0.5, 0.7]
    for threshold in thresholds:
        y_pred_threshold = (y_pred_proba >= threshold).astype(int)
        cm = confusion_matrix(y, y_pred_threshold)
        f1 = f1_score(y, y_pred_threshold)
        wandb.log({f"{dataset_name}_confusion_matrix_threshold_{threshold}": cm, f"{dataset_name}_f1_score_threshold_{threshold}": f1})
        print(f"Confusion Matrix at threshold {threshold} - {dataset_name}:\n{cm}")
        print(f"F1 Score at threshold {threshold} - {dataset_name}: {f1}")
        plt.figure()
        sns.heatmap(cm, annot=True, fmt='d', cmap='Blues')
        plt.title(f'Confusion Matrix at threshold {threshold} - {dataset_name}')
        plt.xlabel('Predicted')
        plt.ylabel('Actual')
        cm_path = os.path.join(plot_output_folder, f'{model_name}_confusion_matrix_{dataset_name}_threshold_{threshold}.png')
        plt.savefig(cm_path)
        plt.close()
        wandb.log({f"{dataset_name}_confusion_matrix_{threshold}": wandb.Image(cm_path)})

if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Train a model for case-control classification.')
    parser.add_argument('--model_name', type=str, required=True, help="Model type: 'logistic_regression', 'random_forest', 'xgboost', 'svm'")
    parser.add_argument('--plot_output_folder', type=str, required=True, help="Folder path to save the plots and model")
    parser.add_argument('--model_output_folder', type=str, required=True, help="Folder path to save the model")
    parser.add_argument('--feature_selection', type=str, required=True, help="Define whether or not feature selection is applied prior to training")
    parser.add_argument('--features_to_bypass_fs', type=str, default='', help="CSV file of features to bypass feature selection")
    parser.add_argument('--features_to_select_fs', type=str, default='', help="CSV file of features to select for feature selection")
    parser.add_argument('--target_variable', type=str, default='', help="The target variable is either prevalent (i.e HCM case/control status) or incident (i.e incident HCM diagnosis)")
    parser.add_argument('--X_train_data', type=str, default='', help="The model data containing training data for X")
    parser.add_argument('--y_train_data', type=str, default='', help="The model data containing training data for y")
    parser.add_argument('--filesuffix', type=str, default='', help="Suffix to append to the model name")
    args = parser.parse_args()

    if args.target_variable in ['allcases','prevalent']: # Initialize wandb with a project name based on the model type
        wandb.init(project=f"polyproteomic_casecontrol_{args.model_name}")
    elif args.target_variable == 'incident':
        wandb.init(project=f"polyproteomic_incident_{args.model_name}")

    # Log configuration parameters
    wandb.config.update({
        "model_type": args.model_name,
        "plot_output_folder": args.plot_output_folder,
        "model_output_folder": args.model_output_folder,
        "feature_selection": args.feature_selection
    })

    #If plotting output folder does not exist, create it
    if not os.path.exists(args.plot_output_folder):
        os.makedirs(args.plot_output_folder)


    print("Importing data...")
    wandb.log({"status": "Importing data"})

    X_train = pd.read_csv(args.X_train_data)
    y_train = pd.read_csv(args.y_train_data)

    #Remove the eid columns from X_train
    if 'eid' in X_train.columns:
        X_train = X_train.drop(columns='eid')

    print(X_train.shape)

    #Print the number of HCM=true and HCM=false in the training and test set
    print(y_train.value_counts())

    print("Data import completed.")
    wandb.log({"status": "Data import completed"})

    #Load the features to bypass and the features to select
    if args.features_to_bypass_fs:
        features_to_bypass_fs = pd.read_csv(args.features_to_bypass_fs, header=None).iloc[:, 0].tolist()
        print('Quantitative features to bypass feature selection')
        print(features_to_bypass_fs)
        wandb.config.update({"features_to_bypass_fs": features_to_bypass_fs})
    else:
        features_to_bypass_fs = []

    if args.features_to_select_fs:
        features_to_select_fs = pd.read_csv(args.features_to_select_fs).loc[:, 'name'].tolist()
        print('Quantitative features of interest to select for')
        print(features_to_select_fs)
        wandb.config.update({"features_to_select_fs": features_to_select_fs})
    else:
        features_to_select_fs = []

    # Get the model and parameter grid based on user input
    model, param_grid = get_model_and_params(args.model_name)

    # Log the parameter grid to wandb config
    wandb.config.update(param_grid)

    if args.filesuffix: #If there is a filesuffix, also add that the model name
        args.model_name = args.model_name + '_' + args.filesuffix

    # Train the model
    trained_model = train_model(X_train, y_train, model, param_grid, args.model_name, args.model_output_folder, args.X_train_data, args.feature_selection, features_to_bypass_fs, features_to_select_fs)

    # Plot metrics for the training set
    plot_metrics(trained_model, X_train, y_train, "Training", args.model_name, args.plot_output_folder)