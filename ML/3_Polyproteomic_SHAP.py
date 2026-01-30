'''
Script to evaluate the SHAP values of the trained model on the entire dataset and output the SHAP summary plot.
Author: Jonathan Chan
Date: 2025-02-28
'''

import os
import joblib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import shap
import argparse
import sklearn
from adjustText import adjust_text
from sklearn.feature_selection import f_classif
from sklearn.utils import resample
import matplotlib.gridspec as gridspec # Needed for custom subplot layout

sklearn.set_config(transform_output="pandas")

def load_model(model_path):
    """Load the saved model from a .pkl file."""
    print(f"Loading model from {model_path}")
    return joblib.load(model_path)

# Filter for specific features by name
def filter_shap_values(shap_values, feature_names, features_to_keep):
    # Get indices of features to keep
    indices = [feature_names.index(feature) for feature in features_to_keep]
    
    # Create a filtered version of the SHAP values
    filtered_values = shap_values.values[:, indices]
    
    # Create a new SHAP values object
    filtered_shap_values = shap.Explanation(
        values=filtered_values,
        base_values=shap_values.base_values,
        data=shap_values.data[:, indices],
        feature_names=[feature_names[i] for i in indices]
    )
    
    return filtered_shap_values

def plot_shap_plots(shap_values, model_name, n_features, output_folder, suffix=''):
    """Plot SHAP plots for the given model."""
    print(f"Plotting SHAP plots for {model_name}")

    #PLOT BAR PLOT
    # Extract feature names directly from shap_values
    feature_names = shap_values.feature_names
    
    # Aggregate mean and standard deviation for each feature
    shap_mean = np.abs(shap_values.values).mean(axis=0)
    shap_std = np.abs(shap_values.values).std(axis=0)
    
    # Sort by mean magnitude
    sorted_idx = np.argsort(shap_mean)[::-1]
    
    # Slice for top n features
    shap_mean_sorted = shap_mean[sorted_idx][:n_features]
    shap_std_sorted = shap_std[sorted_idx][:n_features]
    feature_names_sorted = np.array(feature_names)[sorted_idx][:n_features]
    
    # Plot bar chart with error bars
    plt.figure(figsize=(10, 6))
    plt.barh(feature_names_sorted, shap_mean_sorted, xerr=shap_std_sorted, color='skyblue')
    plt.gca().invert_yaxis()
    plt.xlabel("Mean |SHAP|")
    plt.title(f'Mean Absolute SHAP Values for each of {n_features} features')
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, f"{model_name}_shap_barplot{suffix}.png"))
    plt.close()

    #PLOT BEESWARM PLOT
    shap.plots.beeswarm(shap_values, max_display=n_features, show=False)
    plt.title(f'SHAP beeswarm plot for {n_features} features')
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, f"{model_name}_shap_beeswarm{suffix}.png"))
    plt.close()

def dependence_shap_plotter(shap_values, model_name, output_folder, top_n=5, color_by=None, suffix='', alpha_majority=0.1, marker_size=10):
    """
    Plot SHAP dependence plots for the top_n features.
    Mitigates overplotting for imbalanced binary color_by arrays by making the majority class transparent.
    Includes a marginal histogram of feature values when color_by is specified.

    Parameters:
    - shap_values: SHAP explanation object (e.g., shap.Explanation).
                   Requires shap_values.data, shap_values.values, shap_values.feature_names.
    - model_name: Name of the model for plot titles and filenames.
    - output_folder: Path to save the plots.
    - top_n: Number of top features to plot.
    - color_by: 1D numpy array of binary values (0 or 1) to color points by.
                Assumes 0 is the majority class ('Controls') and 1 is the minority class ('Cases').
                If None, uses default SHAP interaction coloring and histogram.
    - suffix: Optional string to append to the output filename.
    - alpha_majority: Alpha transparency for the 'Controls' class (label 0). Default 0.1.
                      Adjust based on imbalance ratio.
    - marker_size: Size of the scatter plot markers. Default 10.
    """
    # Ensure output folder exists
    os.makedirs(output_folder, exist_ok=True)

    # --- Data Extraction and Validation ---
    if isinstance(shap_values, shap.Explanation):
        shap_vals_arr = shap_values.values
        feature_names = shap_values.feature_names
        base_data = shap_values.data # This is the original feature data
    else:
        raise TypeError("Expected shap_values to be a shap.Explanation object containing .values, .data, and .feature_names")

    if feature_names is None:
         raise ValueError("shap_values object must have feature_names attribute.")
    if base_data is None:
         raise ValueError("shap_values object must have data attribute (original features).")

    # --- Feature Importance and Selection ---
    feature_importance = np.abs(shap_vals_arr).mean(0)
    # Ensure top_n does not exceed the number of features
    actual_top_n = min(top_n, len(feature_names))
    if actual_top_n < top_n:
        print(f"Warning: Requested top_n={top_n}, but only {len(feature_names)} features available. Plotting top {actual_top_n}.")
    top_indices = np.argsort(-feature_importance)[:actual_top_n] # Get indices of top features

    # --- Plotting Setup ---
    ncols = 3
    nrows = (actual_top_n + ncols - 1) // ncols # Calculate rows needed
    fig = plt.figure(figsize=(18, 6 * nrows)) # Create figure directly
    # Define the outer GridSpec for the entire figure
    outer_grid = gridspec.GridSpec(nrows, ncols, figure=fig)

    plot_index = 0 # Keep track of which subplot position we are using

    # --- Plotting Loop ---
    for i, feature_idx in enumerate(top_indices):
        # Get the SubplotSpec for the current plot's cell in the outer grid
        subplot_spec = outer_grid[plot_index]

        feature_name = feature_names[feature_idx]

        # Get SHAP values and feature values for the current feature
        current_shap_values = shap_vals_arr[:, feature_idx]

        # Extract corresponding feature values from base_data
        if isinstance(base_data, pd.DataFrame):
             if feature_name not in base_data.columns:
                 raise ValueError(f"Feature '{feature_name}' not found in base_data columns.")
             current_feature_values = base_data[feature_name].values
        elif isinstance(base_data, np.ndarray):
             if feature_idx >= base_data.shape[1]:
                  raise IndexError(f"feature_idx {feature_idx} out of bounds for base_data with shape {base_data.shape}")
             current_feature_values = base_data[:, feature_idx]
        else:
             raise TypeError(f"Unsupported type for base_data: {type(base_data)}. Expected pandas DataFrame or NumPy array.")

        # --- Conditional Plotting ---
        if color_by is not None:
            # --- Manual Plotting with Histogram (color_by provided) ---

            # Create a nested GridSpec *within* the current subplot cell
            # 2 rows, 1 column. Scatter plot taller than histogram.
            nested_gs = gridspec.GridSpecFromSubplotSpec(
                2, 1, subplot_spec=subplot_spec, height_ratios=[4, 1], hspace=0.05
            )

            # Add subplots using the nested GridSpec indices
            ax_scatter = fig.add_subplot(nested_gs[0]) # Axes for scatter plot (top row)
            ax_hist = fig.add_subplot(nested_gs[1], sharex=ax_scatter) # Axes for histogram (bottom row)

            # Validate color_by array
            if not isinstance(color_by, np.ndarray):
                try:
                    color_by = np.array(color_by)
                except Exception as e:
                     raise TypeError(f"color_by must be convertible to a numpy array. Error: {e}")

            if color_by.ndim != 1 or len(color_by) != len(current_shap_values):
                raise ValueError(f"color_by must be a 1D array with the same length as the number of samples ({len(current_shap_values)}). Found shape {color_by.shape}")

            unique_labels = np.unique(color_by)
            if not np.all(np.isin(unique_labels, [0, 1])):
               print(f"Warning: color_by contains labels other than 0 and 1 ({unique_labels}). Plotting may not be as intended.")

            # Identify indices for Controls (0) and Cases (1)
            idx_controls = np.where(color_by == 0)[0]
            idx_cases = np.where(color_by == 1)[0]

            # Plot Controls points first with transparency
            ax_scatter.scatter(current_feature_values[idx_controls],
                       current_shap_values[idx_controls],
                       color='blue', # Or choose another color
                       alpha=alpha_majority,
                       s=marker_size,
                       label='Controls', # Updated label
                       rasterized=True)

            # Plot Cases points second, fully opaque
            ax_scatter.scatter(current_feature_values[idx_cases],
                       current_shap_values[idx_cases],
                       color='red', # Or choose another color
                       alpha=1.0,
                       s=marker_size,
                       label='Cases', # Updated label
                       rasterized=True)

            ax_scatter.legend() # Add legend to distinguish classes

            # --- Scatter Plot Styling ---
            ax_scatter.set_title(f'Feature: {feature_name}', fontsize=12)
            ax_scatter.set_ylabel(f'SHAP value for {feature_name}')
            ax_scatter.axhline(y=0, color='black', linestyle='--', linewidth=0.5)
            # Hide x-axis labels and ticks on the scatter plot
            plt.setp(ax_scatter.get_xticklabels(), visible=False)
            ax_scatter.tick_params(axis='x', which='both', length=0) # Hide x-ticks


            # --- Histogram Plotting ---
            # Plot histogram of the feature values on the bottom axes
            ax_hist.hist(current_feature_values, bins=50, color='lightgray', density=True, align='mid')

            # --- Histogram Styling ---
            ax_hist.set_xlabel(f'{feature_name} Value')
            # Make histogram less visually intrusive
            ax_hist.set_yticks([])
            ax_hist.set_yticklabels([])
            ax_hist.spines['top'].set_visible(False)
            ax_hist.spines['right'].set_visible(False)
            ax_hist.spines['left'].set_visible(False)
            ax_hist.tick_params(axis='x', direction='in') # Ticks inside

            # hspace is now set in GridSpecFromSubplotSpec

        else:
            # --- Default SHAP Plotting (color_by is None) ---
            # Use the standard shap plot which includes interaction coloring and histogram
            # Create a single axes occupying the whole subplot cell
            ax = fig.add_subplot(subplot_spec)

            # Construct the Explanation object slice needed by shap.plots.scatter
            shap_values_slice = shap_values[:, feature_idx]

            shap.plots.scatter(
                shap_values_slice, # Pass the explanation slice for this feature
                color=None, # Let SHAP choose interaction feature
                show=False,
                ax=ax
            )
            # Standard plot already includes title, labels, line, and histogram

        plot_index += 1 # Move to the next subplot position

    # --- Final Figure Adjustments ---
    # Add an overall title to the figure
    fig.suptitle(f'SHAP Dependence Plots for Top {actual_top_n} Features ({model_name})', fontsize=16, y=1.0)
    # Adjust layout - use tight_layout first, then potentially refine with subplots_adjust
    # Use wspace and hspace in the outer_grid definition if needed, or subplots_adjust
    outer_grid.tight_layout(fig, rect=[0, 0.03, 1, 0.97]) # Apply tight_layout to the grid

    # Save the figure
    output_filename = os.path.join(output_folder, f'{model_name}_top{actual_top_n}_dependence_plots{suffix}.png')
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"Saved dependence plots to {output_filename}")
    plt.close(fig) # Close the figure to free memory


def mean_shap_to_df(shap_values):
    # Calculate mean absolute SHAP value for each feature
    mean_abs_shap = np.abs(shap_values.values).mean(axis=0)
    
    # Create DataFrame with feature names and their importance
    feature_importance = pd.DataFrame({
        'Feature': shap_values.feature_names,
        'mean_abs_shap_value': mean_abs_shap
    })
    
    # Sort by importance (highest to lowest)
    feature_importance = feature_importance.sort_values('mean_abs_shap_value', ascending=False)
    
    return feature_importance

from typing import List, Optional, Union, Any
def save_precomputed_shap_values(
    shap_values_data: Any, # Pre-computed SHAP values (Explanation obj, np.ndarray, or list)
    feature_names: List[str], # List of feature names (now mandatory)
    output_path: str, # Directory path to save the output files
    model_name: str, # Descriptive name for the model, used in filenames
    suffix: Optional[str] = None, # Optional suffix string for filenames
    base_value: Optional[Union[float, np.ndarray, list]] = None, # Provide if not in shap_values_data or to override
    instance_index: Optional[pd.Index] = None, # Optional: provide index for rows
    positive_class_index: Optional[int] = 1 # For multi-output classification data
):
    """
    Saves pre-computed SHAP values (and optionally base value) to TSV/TXT files.

    Accepts SHAP values as a SHAP Explanation object, a NumPy array, or a list
    of NumPy arrays (for multi-output).

    Args:
        shap_values_data: The pre-computed SHAP values. Can be:
                          - A SHAP Explanation object (`explanation = explainer(X)`).
                          - A NumPy array (n_instances, n_features) for single output.
                          - A list of NumPy arrays [(n_instances, n_features)] for multi-output.
        feature_names: A list of strings representing the feature names. The order
                       must match the feature dimension of the SHAP values. Mandatory.
        output_path: Directory path where the output files will be saved.
                     The directory will be created if it doesn't exist.
        model_name: Descriptive name for the model, used as the base for filenames.
        suffix: An optional suffix string to append to the filenames before the
                extension (e.g., '_test_set'). An underscore is added automatically.
        base_value: The base value (or list/array for multi-output) corresponding
                    to the SHAP values. If `shap_values_data` is an Explanation
                    object, its `base_values` are used by default unless this argument
                    is provided (override). If `shap_values_data` is array/list,
                    this argument is the only way to save the base value.
        instance_index: Optional Pandas Index object to use for the rows in the
                        output TSV file. If None, no index is saved.
        positive_class_index: Applicable if `shap_values_data` represents multi-output
                               classification values (e.g., a list of arrays or a 3D array).
                               Specifies the class index to extract and save.
                               If None and multi-output is detected, raises ValueError.
                               Ignored for single-output 2D arrays.

    Raises:
        ValueError: If dimensions mismatch, index mismatch, or required info is missing.
        TypeError: If `shap_values_data` is an unexpected type.
        RuntimeError: If SHAP values array cannot be processed.
    """
    print(f"Saving pre-computed SHAP values for: {model_name}{'_' + suffix if suffix else ''}")

    # --- 1. Interpret Input & Extract Raw SHAP Values ---
    shap_values_raw: Any
    _base_value = base_value # Start with the provided base_value
    _feature_names = feature_names # Use the mandatory feature_names

    # Check if input is a SHAP Explanation object
    is_explanation_obj = hasattr(shap_values_data, 'values') and \
                         hasattr(shap_values_data, 'base_values')
                         # We don't strictly need .feature_names here as it's mandatory arg

    if is_explanation_obj:
        print("Input detected as SHAP Explanation object.")
        shap_values_raw = shap_values_data.values
        if _base_value is None: # If user didn't override, use the object's base value
            _base_value = shap_values_data.base_values
            print("Using base_values found in Explanation object.")
        else:
            print("Using provided base_value (overriding Explanation object's base_values).")

        # Sanity check feature names length if possible
        if hasattr(shap_values_data, 'feature_names') and shap_values_data.feature_names is not None:
            if len(shap_values_data.feature_names) != len(_feature_names):
                print(f"Warning: Provided feature_names count ({len(_feature_names)}) differs from Explanation object's feature_names count ({len(shap_values_data.feature_names)}). Using provided list.")
            # Potentially compare names if needed, but trust the mandatory argument for now.

    elif isinstance(shap_values_data, (np.ndarray, list)):
        print("Input detected as NumPy array or list.")
        shap_values_raw = shap_values_data
        if _base_value is None:
            print("Warning: No base_value provided for raw SHAP input. Base value file will not be saved.")
    else:
        raise TypeError(f"Unsupported type for shap_values_data: {type(shap_values_data)}. Expected Explanation object, np.ndarray, or list.")

    # --- 2. Handle Multi-Output Structure ---
    shap_values_array: np.ndarray
    n_instances = None # Track number of instances

    # Determine if the raw values structure indicates multiple outputs
    is_multi_output = isinstance(shap_values_raw, list) or \
                      (isinstance(shap_values_raw, np.ndarray) and shap_values_raw.ndim > 2)

    if is_multi_output:
         # Try to determine number of outputs
         if isinstance(shap_values_raw, list):
             if not shap_values_raw: raise ValueError("Input SHAP values list is empty.")
             num_outputs = len(shap_values_raw)
         elif isinstance(shap_values_raw, np.ndarray):
             # Guess output dim: often last for SHAP (n_inst, n_feat, n_out) or first (n_out, n_inst, n_feat)
             # This is ambiguous without more context, we rely on positive_class_index slicing logic
             # Let's just proceed and let slicing + validation handle it.
             num_outputs = 'unknown (3D array)' # Placeholder

         print(f"Multi-output SHAP values detected (structure suggests {num_outputs} outputs).")

         if positive_class_index is None:
             raise ValueError(f"Multi-output detected, but 'positive_class_index' is None. Please specify which class index to save.")
         if not isinstance(positive_class_index, int) or positive_class_index < 0:
              raise ValueError(f"'positive_class_index' ({positive_class_index}) must be a non-negative integer.")

         print(f"Selecting SHAP values for class index: {positive_class_index}")
         try:
             if isinstance(shap_values_raw, list):
                 if positive_class_index >= len(shap_values_raw):
                      raise IndexError(f"positive_class_index {positive_class_index} out of bounds for list of length {len(shap_values_raw)}.")
                 shap_values_array = shap_values_raw[positive_class_index]
             else: # Assuming 3D numpy array
                # Try common slicings, prioritize shape that matches feature count
                if shap_values_raw.ndim == 3:
                    if shap_values_raw.shape[1] == len(_feature_names): # Assume (n_instances, n_features, n_outputs)
                        if positive_class_index >= shap_values_raw.shape[2]: raise IndexError("Index out of bounds for last dimension.")
                        shap_values_array = shap_values_raw[:, :, positive_class_index]
                    elif shap_values_raw.shape[2] == len(_feature_names): # Assume (n_instances, n_outputs, n_features)? Less common
                        if positive_class_index >= shap_values_raw.shape[1]: raise IndexError("Index out of bounds for middle dimension.")
                        shap_values_array = shap_values_raw[:, positive_class_index, :]
                    elif shap_values_raw.shape[0] == len(_feature_names): # Assume (n_features, n_instances, n_outputs)? Very unlikely
                        raise ValueError("Cannot reliably interpret 3D array shape {shap_values_raw.shape} with features as first dimension.")
                    else: # Maybe (n_outputs, n_instances, n_features)?
                        if positive_class_index >= shap_values_raw.shape[0]: raise IndexError("Index out of bounds for first dimension.")
                        shap_values_array = shap_values_raw[positive_class_index, :, :]

                else:
                     raise ValueError(f"Cannot handle multi-output NumPy array with shape {shap_values_raw.shape}")

         except IndexError as e:
              raise ValueError(f"Could not slice multi-output SHAP data at index {positive_class_index}: {e}")

         # Try to select the corresponding base value if _base_value is also multi-output
         if _base_value is not None and hasattr(_base_value, '__len__') and not isinstance(_base_value, str):
             try:
                 # Check if index is valid for base value length
                 if positive_class_index < len(_base_value):
                     _base_value = _base_value[positive_class_index]
                     print(f"Selected corresponding base value for class index {positive_class_index}.")
                 else:
                      print(f"Warning: positive_class_index {positive_class_index} is out of bounds for base_value length {len(_base_value)}. Keeping original base value.")
             except (IndexError, TypeError): # Handle non-sliceable or index error
                 print(f"Warning: Could not select base value for index {positive_class_index}. Keeping original base value.")
         elif _base_value is not None:
             print("Warning: SHAP values seem multi-output, but base value does not. Base value might be average.")

    else:
        # Assuming single output (should be a 2D NumPy array)
        print("Assuming single-output SHAP values.")
        if not isinstance(shap_values_raw, np.ndarray):
             raise TypeError(f"Expected single-output SHAP values to be ndarray, but got {type(shap_values_raw)}")
        if shap_values_raw.ndim != 2:
            raise ValueError(f"Expected single-output SHAP values to be 2D (instances, features), but got shape {shap_values_raw.shape}.")
        shap_values_array = shap_values_raw
        # Base value might still be multi-output (e.g., one per class)
        if _base_value is not None and hasattr(_base_value, '__len__') and not isinstance(_base_value, str) and len(_base_value) > 1:
             print(f"Note: Base value has multiple ({len(_base_value)}) entries but SHAP values are single-output. Saving all base values.")

    # --- 3. Final Validation ---
    if not isinstance(shap_values_array, np.ndarray) or shap_values_array.ndim != 2:
         raise RuntimeError(f"Failed to extract a 2D NumPy array for SHAP values. Got shape {getattr(shap_values_array, 'shape', 'N/A')}.")

    n_instances = shap_values_array.shape[0]
    n_features = shap_values_array.shape[1]

    if n_features != len(_feature_names):
         raise ValueError(f"Dimension mismatch: Final SHAP values feature count ({n_features}) != provided feature_names count ({len(_feature_names)}).")
    if instance_index is not None and len(instance_index) != n_instances:
         raise ValueError(f"Dimension mismatch: instance_index length ({len(instance_index)}) != SHAP values instance count ({n_instances}).")

    # --- 4. Create DataFrame ---
    print(f"Creating DataFrame ({n_instances} instances, {n_features} features)...")
    shap_df = pd.DataFrame(shap_values_array, columns=_feature_names, index=instance_index)

    # --- 5. Prepare Filenames & Save ---
    os.makedirs(output_path, exist_ok=True) # Ensure directory exists

    file_suffix_str = f"_{suffix}" if suffix else ""
    shap_values_filename = os.path.join(output_path, f"{model_name}{file_suffix_str}_shap_values.tsv")
    base_value_filename = os.path.join(output_path, f"{model_name}{file_suffix_str}_shap_base_value.txt")

    # Save SHAP values DataFrame
    print(f"Saving SHAP values to: {shap_values_filename}")
    shap_df.to_csv(shap_values_filename, sep='\t', index=(instance_index is not None), float_format='%.8g') # Save index only if provided

    # Save Base value
    if _base_value is not None:
        print(f"Saving base value(s) to: {base_value_filename}")
        try:
            with open(base_value_filename, 'w') as f:
                # Handle array/list/single value cleanly
                if hasattr(_base_value, '__iter__') and not isinstance(_base_value, str):
                    f.write('\t'.join(map(lambda x: f"{x:.8g}", _base_value))) # Format numbers
                else:
                    f.write(f"{_base_value:.8g}") # Format single number
        except Exception as e:
            print(f"Warning: Could not save base value to {base_value_filename}: {e}")
    else:
        print("Base value was not available or not provided, skipping save.")

    print(f"SHAP saving process complete for {model_name}{file_suffix_str}.")

def plot_shap_vs_Fratio(mean_shap_df, fratio_df, model_name, output_folder, top_n_label=5):
    """Plot mean absolute SHAP-values vs F-ratio for the given model. It also labels the top 5 features for both SHAP and F ratio"""
    print(f"Plotting mean absolute SHAP-values vs F-ratio for {model_name}")

    # Merge the feature importances and F-ratio dataframes
    merged_df = mean_shap_df.merge(fratio_df, on='Feature', how='inner')

    #Add the label column with True if in top 5 features for either SHAP or Fratio
    top5_fratio = fratio_df.sort_values('F-ratio',ascending=False).iloc[:top_n_label+1,:]['Feature'].values
    top5_shap = mean_shap_df.sort_values('mean_abs_shap_value',ascending=False).iloc[:top_n_label+1,:]['Feature'].values

    merged_df['label'] = merged_df['Feature'].apply(
        lambda x: x if (x in top5_fratio or x in top5_shap) else '')
    
    plt.figure(figsize=(12, 9))
    sns.scatterplot(x='mean_abs_shap_value', y='F-ratio', data=merged_df)

    # Collect texts
    texts = []
    for i, row in merged_df.iterrows():
        texts.append(plt.text(row['mean_abs_shap_value'], row['F-ratio'], row['label']))
    
    # Adjust label positions to minimize overlaps
    adjust_text(texts)
    
    plt.title(f'Mean Absolute SHAP value vs F-ratio - {model_name}')
    plt.xlabel('Mean Absolute SHAP value')
    plt.ylabel('F-ratio')
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, f'{model_name}_shap_vs_fratio.png'))
    # plt.show()
    plt.close()
    print(f"Feature importance vs F-ratio plot saved for {model_name}")

##############################################
# New functions for SHAP interaction analysis #
##############################################
def filter_shap_interaction_values(interaction_values, feature_names, features_to_keep):
    """
    Filter the SHAP interaction values to only include the specified features.
    interaction_values: numpy array with shape (n_samples, n_features, n_features)
    feature_names: list of feature names corresponding to the second and third axes
    features_to_keep: list of feature names to retain
    """
    indices = [feature_names.index(feature) for feature in features_to_keep]
    # Filter along both feature dimensions
    filtered_interaction = interaction_values[:, indices, :][:, :, indices]
    return filtered_interaction

def plot_shap_interaction_heatmap(mean_interaction, feature_names, model_name, output_folder, suffix=''):
    """
    Plot a heatmap of the mean absolute SHAP interaction values.
    mean_interaction: 2D numpy array of shape (n_features, n_features)
    feature_names: list of feature names corresponding to the matrix dimensions
    """
    plt.figure(figsize=(12, 10))
    sns.heatmap(mean_interaction, xticklabels=feature_names, yticklabels=feature_names,
                cmap='viridis', annot=True, fmt=".3f")
    plt.title(f"Mean Absolute SHAP Interaction Values - {model_name}")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, f"{model_name}_shap_interaction_heatmap{suffix}.png"))
    plt.close()

##############################################
# New functions for SHAP Bootstrapping #
##############################################

def bootstrap_shap_from_X(X, model, n_iterations=1000, alpha=0.05):
    """
    Perform bootstrapping by resampling X_train and computing SHAP values for each bootstrap sample.
    Returns mean, lower bound, and upper bound for each feature.
    """
    n_samples = X.shape[0]
    n_features = X.shape[1]
    bootstrap_means = np.zeros((n_iterations, n_features))
    
    explainer = shap.Explainer(model.named_steps['classifier'])
    
    for i in range(n_iterations):
        X_resampled = resample(X, replace=True, n_samples=n_samples, random_state=i)
        shap_explainer_obj = explainer(X_resampled)
        shap_values = shap_explainer_obj.values
        bootstrap_means[i, :] = np.abs(shap_values).mean(axis=0)
    
    mean_vals = np.mean(bootstrap_means, axis=0)
    lower_bound = np.percentile(bootstrap_means, 100 * (alpha / 2), axis=0)
    upper_bound = np.percentile(bootstrap_means, 100 * (1 - alpha / 2), axis=0)

    # Extract feature names directly from shap_values
    feature_names = shap_explainer_obj.feature_names
    
    return mean_vals, lower_bound, upper_bound, feature_names

def plot_bootstrap_shap(mean_vals, lower_bound, upper_bound, feature_names, n_features, model_name, output_folder, suffix=''):
    """
    Plot the bootstrapped SHAP values with confidence intervals.
    """
    # Sort by mean magnitude
    sorted_idx = np.argsort(mean_vals)[::-1]

    #Filter mean_vals, lower_bound, upper_bound, feature_names to n_features
    mean_vals = mean_vals[sorted_idx][:n_features]
    lower_bound = lower_bound[sorted_idx][:n_features]
    upper_bound = upper_bound[sorted_idx][:n_features]
    feature_names = np.array(feature_names)[sorted_idx][:n_features]

    plt.figure(figsize=(10, 6))
    plt.barh(feature_names, mean_vals, xerr=[mean_vals - lower_bound, upper_bound - mean_vals], color='skyblue')
    plt.gca().invert_yaxis()
    plt.xlabel("Mean |SHAP|")
    plt.title(f'Bootstrapped Mean Absolute SHAP Values for {model_name}')
    plt.tight_layout()
    plt.savefig(os.path.join(output_folder, f"{model_name}_bootstrapped_shap{suffix}.png"))
    # plt.show()
    plt.close()

    print(f"Bootstrapped SHAP values plot saved for {model_name}")

#############################################################################################
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Evaluate a saved model on a test set and output performance plots.')
    parser.add_argument('--model_pkl_file', type=str, required=True, help='Path to the saved model .pkl file')
    parser.add_argument('--X_train_preprocessed_path', type=str, required=True, help='Path to the X_train_preprocessed file')
    parser.add_argument('--y_train_path',type=str, required=True, help='Path to the y_train file')
    parser.add_argument('--plot_output_path', type=str, required=True, help='Path to save the output plots')
    parser.add_argument('--model_name', type=str, required=True, help='Name of the model for plot titles and filenames')
    parser.add_argument('--pp_names_file', type=str, required=True, help="CSV file containing the plasma protein names")
    parser.add_argument('--bootstrap_or_not', type=str, default='no', help='Whether to compute bootstrapped SHAP values or not. Default is "no".')
    args = parser.parse_args()

    ##############################################
    # Preparation and computation of SHAP values  #
    ##############################################

    # Create the output folder if it does not exist
    if not os.path.exists(args.plot_output_path):
        os.makedirs(args.plot_output_path)
        print(f"Created output folder: {args.plot_output_path}")

    # Load the data
    print(f"Loading data from {args.X_train_preprocessed_path}")
    X = pd.read_csv(args.X_train_preprocessed_path)
    y = pd.read_csv(args.y_train_path).values.ravel()
    print(f"Preprocessed training data loaded with shape: {X.shape}")

    pp_names = pd.read_csv(args.pp_names_file).iloc[:, 0].values

    model = load_model(args.model_pkl_file)
    explainer = shap.Explainer(model.named_steps['classifier'])

    # Separate cases and controls
    cases = X[y == 1]
    controls = X[y == 0]

    # Compute SHAP values for cases and controls
    print("Computing SHAP values for cases...")
    shap_values_cases = explainer(cases)
    print("Computing SHAP values for controls...")
    shap_values_controls = explainer(controls)

    ##############################################
    # Compute and plot basic SHAP plots comparing features  #
    ##############################################

    # Plot SHAP summary plots for both cases and controls combined
    shap_values_combined = explainer(X)

    if len(X.columns) <= 50:
        plot_shap_plots(shap_values_combined, args.model_name, len(X.columns), args.plot_output_path, '_combined')
    plot_shap_plots(shap_values_combined, args.model_name, 10, args.plot_output_path, '_combined_top10')

    # Plot SHAP summary plots for cases only
    if len(cases.columns) <= 50:
        plot_shap_plots(shap_values_cases, args.model_name, len(cases.columns), args.plot_output_path, '_cases')
    plot_shap_plots(shap_values_cases, args.model_name, 10, args.plot_output_path, '_cases_top10')

    # Plot SHAP summary plots for controls only
    if len(controls.columns) <= 50:
        plot_shap_plots(shap_values_controls, args.model_name, len(controls.columns), args.plot_output_path, '_controls')
    plot_shap_plots(shap_values_controls, args.model_name, 10, args.plot_output_path, '_controls_top10')

    # Filter SHAP values for plasma proteins
    pp_names = pp_names[np.isin(pp_names, X.columns.values)]
    print(f'The filtered features include {pp_names}')
    shap_values_cases_filtered = filter_shap_values(shap_values_cases, X.columns.values.tolist(), pp_names)
    shap_values_controls_filtered = filter_shap_values(shap_values_controls, X.columns.values.tolist(), pp_names)
    shap_values_combined_filtered = filter_shap_values(shap_values_combined, X.columns.values.tolist(), pp_names)

    #Save the filtered SHAP value objects
    save_precomputed_shap_values(shap_values_cases_filtered, pp_names, args.plot_output_path, args.model_name, suffix='_cases_ppfiltered')
    save_precomputed_shap_values(shap_values_controls_filtered, pp_names, args.plot_output_path, args.model_name, suffix='_controls_ppfiltered')
    save_precomputed_shap_values(shap_values_combined_filtered, pp_names, args.plot_output_path, args.model_name, suffix='_combined_ppfiltered')

    #Also replot all the SHAP summary plots for the filtered features
    plot_shap_plots(shap_values_cases_filtered, args.model_name, len(pp_names), args.plot_output_path, '_cases_ppfiltered')
    plot_shap_plots(shap_values_controls_filtered, args.model_name, len(pp_names), args.plot_output_path, '_controls_ppfiltered')
    plot_shap_plots(shap_values_combined_filtered, args.model_name, len(pp_names), args.plot_output_path, '_combined_ppfiltered')

    # Plot SHAP dependence plots for filtered features
    dependence_shap_plotter(shap_values_cases_filtered, args.model_name, args.plot_output_path, top_n=9,suffix='_cases_ppfiltered')
    dependence_shap_plotter(shap_values_controls_filtered, args.model_name, args.plot_output_path, top_n=9,suffix='_controls_ppfiltered')
    dependence_shap_plotter(shap_values_combined_filtered, args.model_name, args.plot_output_path, top_n=9, color_by=y,suffix='_combined_ppfiltered')

    ##############################################
    # Compute and plot SHAP values against the F-ratio   #
    ##############################################

    # Compute F-ratio for plasma proteins
    fvalues = f_classif(X.loc[:, pp_names], y)[0]
    fratio_df = pd.DataFrame({'Feature': pp_names, 'F-ratio': fvalues})

    # Compute mean SHAP values for cases and controls
    mean_shap_cases_df = mean_shap_to_df(shap_values_cases_filtered)
    mean_shap_controls_df = mean_shap_to_df(shap_values_controls_filtered)
    mean_shap_combined_df = mean_shap_to_df(shap_values_combined)

    # Save mean SHAP values
    mean_shap_cases_df.to_csv(os.path.join(args.plot_output_path, f'{args.model_name}_mean_abs_shap_values_cases.csv'), index=False)
    mean_shap_controls_df.to_csv(os.path.join(args.plot_output_path, f'{args.model_name}_mean_abs_shap_values_controls.csv'), index=False)
    mean_shap_combined_df.to_csv(os.path.join(args.plot_output_path, f'{args.model_name}_mean_abs_shap_values_combined.csv'), index=False)

    # Plot SHAP vs F-ratio for cases and controls
    plot_shap_vs_Fratio(mean_shap_cases_df, fratio_df, f"{args.model_name}_cases", args.plot_output_path)
    plot_shap_vs_Fratio(mean_shap_controls_df, fratio_df, f"{args.model_name}_controls", args.plot_output_path)
    plot_shap_vs_Fratio(mean_shap_combined_df, fratio_df, f"{args.model_name}_combined", args.plot_output_path)

    ##############################################
    # Compute and plot SHAP interaction values   #
    ##############################################

    print("Computing SHAP interaction values for cases...")
    shap_interaction_values_cases = explainer.shap_interaction_values(cases)
    print("Computing SHAP interaction values for controls...")
    shap_interaction_values_controls = explainer.shap_interaction_values(controls)

    shap_interaction_values_combined = explainer.shap_interaction_values(X)

    # Filter interaction values for plasma proteins
    shap_interaction_values_cases_filtered = filter_shap_interaction_values(
        shap_interaction_values_cases, X.columns.values.tolist(), pp_names)
    shap_interaction_values_controls_filtered = filter_shap_interaction_values(
        shap_interaction_values_controls, X.columns.values.tolist(), pp_names)
    shap_interaction_values_combined_filtered = filter_shap_interaction_values(
        shap_interaction_values_combined, X.columns.values.tolist(), pp_names 
    )

    save_precomputed_shap_values(shap_interaction_values_cases_filtered, pp_names, args.plot_output_path, args.model_name, suffix='interaction_cases_ppfiltered')
    save_precomputed_shap_values(shap_interaction_values_controls_filtered, pp_names, args.plot_output_path, args.model_name, suffix='interaction_controls_ppfiltered')
    save_precomputed_shap_values(shap_interaction_values_combined_filtered, pp_names, args.plot_output_path, args.model_name, suffix='interaction_combined_ppfiltered')

    # Compute mean absolute interaction values
    mean_interaction_cases = np.abs(shap_interaction_values_cases_filtered).mean(axis=0)
    mean_interaction_controls = np.abs(shap_interaction_values_controls_filtered).mean(axis=0)
    mean_interaction_combined = np.abs(shap_interaction_values_combined_filtered).mean(axis=0)

    # Plot heatmaps for interaction effects
    plot_shap_interaction_heatmap(mean_interaction_cases, pp_names, f"{args.model_name}_cases", args.plot_output_path, suffix='_ppfiltered')
    plot_shap_interaction_heatmap(mean_interaction_controls, pp_names, f"{args.model_name}_controls", args.plot_output_path, suffix='_ppfiltered')
    plot_shap_interaction_heatmap(mean_interaction_combined, pp_names, f"{args.model_name}_combined", args.plot_output_path, suffix='_ppfiltered')

    ##############################################
    # Compute and plot SHAP bootstrap distributions for 95% confidence interval estimation of per-feature mean |SHAP|   #
    ##############################################

    if args.bootstrap_or_not == 'yes':
        print("Computing bootstrapped SHAP values for cases...")
        mean_vals_cases, lower_bound_cases, upper_bound_cases, feature_names_cases = bootstrap_shap_from_X(cases, model, n_iterations=1000, alpha=0.05)
        bootstrap_shap_cases_df = pd.DataFrame({'Feature': feature_names_cases, 'Mean SHAP': mean_vals_cases, 'Lower Bound': lower_bound_cases, 'Upper Bound': upper_bound_cases})
        bootstrap_shap_cases_df.to_csv(os.path.join(args.plot_output_path, f'{args.model_name}_bootstrapped_shap_values_cases.csv'), index=False)
        plot_bootstrap_shap(mean_vals_cases, lower_bound_cases, upper_bound_cases, feature_names_cases, len(feature_names_cases), f"{args.model_name}_cases", args.plot_output_path)
        feature_names_cases_filtered = np.array(feature_names_cases)[np.isin(feature_names_cases, pp_names)]
        plot_bootstrap_shap(mean_vals_cases[np.isin(feature_names_cases, pp_names)], lower_bound_cases[np.isin(feature_names_cases, pp_names)], upper_bound_cases[np.isin(feature_names_cases, pp_names)], 
                            feature_names_cases_filtered, len(feature_names_cases_filtered), f"{args.model_name}_cases", args.plot_output_path, suffix='ppfiltered')

        print("Computing bootstrapped SHAP values for controls...")
        mean_vals_controls, lower_bound_controls, upper_bound_controls, feature_names_controls = bootstrap_shap_from_X(controls, model, n_iterations=1000, alpha=0.05)
        bootstrap_shap_controls_df = pd.DataFrame({'Feature': feature_names_controls, 'Mean SHAP': mean_vals_controls, 'Lower Bound': lower_bound_controls, 'Upper Bound': upper_bound_controls})
        bootstrap_shap_controls_df.to_csv(os.path.join(args.plot_output_path, f'{args.model_name}_bootstrapped_shap_values_controls.csv'), index=False)
        plot_bootstrap_shap(mean_vals_controls, lower_bound_controls, upper_bound_controls, feature_names_controls, len(feature_names_controls), f"{args.model_name}_controls", args.plot_output_path)
        feature_names_controls_filtered= np.array(feature_names_controls)[np.isin(feature_names_controls, pp_names)]
        plot_bootstrap_shap(mean_vals_controls[np.isin(feature_names_controls, pp_names)], lower_bound_controls[np.isin(feature_names_controls, pp_names)], upper_bound_controls[np.isin(feature_names_controls, pp_names)], 
                    feature_names_controls_filtered, len(feature_names_controls_filtered), f"{args.model_name}_controls", args.plot_output_path, suffix='ppfiltered')


        print('Computing bootstrapped SHAP values for combined...')
        mean_vals_combined, lower_bound_combined, upper_bound_combined, feature_names_combined = bootstrap_shap_from_X(X, model, n_iterations=1000, alpha=0.05)
        # Save bootstrapped SHAP values
        bootstrap_shap_combined_df = pd.DataFrame({'Feature': feature_names_combined, 'Mean SHAP': mean_vals_combined, 'Lower Bound': lower_bound_combined, 'Upper Bound': upper_bound_combined})
        bootstrap_shap_combined_df.to_csv(os.path.join(args.plot_output_path, f'{args.model_name}_bootstrapped_shap_values_combined.csv'), index=False)
        # Plot bootstrapped SHAP values
        plot_bootstrap_shap(mean_vals_combined, lower_bound_combined, lower_bound_combined, feature_names_combined, len(feature_names_combined), f"{args.model_name}_combined", args.plot_output_path)
        #For each of the cases, controls and combined, filter the bootstrapped SHAP values to only include the plasma proteins
        feature_names_combined_filtered = np.array(feature_names_combined)[np.isin(feature_names_combined, pp_names)]
        # Plot the bootstrapped SHAP values with confidence intervals
        plot_bootstrap_shap(mean_vals_combined[np.isin(feature_names_combined, pp_names)], lower_bound_combined[np.isin(feature_names_combined, pp_names)], upper_bound_combined[np.isin(feature_names_combined, pp_names)],
                            feature_names_combined_filtered, len(feature_names_combined_filtered), f"{args.model_name}_combined", args.plot_output_path, suffix='ppfiltered')

