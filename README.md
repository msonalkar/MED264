Thymoma Classification Analysis

This repository contains code for analyzing and classifying thymoma cancer subtypes using machine learning approaches on DNA methylation data. The analysis implements multiple classification models to predict thymoma subtypes based on beta values from methylation arrays.

Overview

The project focuses on developing supervised learning models for thymoma classification using the following approaches:

- Support Vector Classification (SVC)
- Extreme Gradient Boosting (XGBoost)
- Multilayer Perceptron Neural Networks (MLP)

Two dimensionality reduction techniques are employed to handle the high-dimensional methylation data:

- Principal Component Analysis (PCA) with 20 and 200 components
- SelectKBest feature selection with 300 features

Requirements

The code requires the following Python packages:

- numpy
- pandas 
- scikit-learn
- xgboost
- tensorflow
- matplotlib
- seaborn
- autograd

Data

The analysis uses two main data files:

- `combinedBetas_batch_corrected.csv` - DNA methylation beta values
- `combinedPheno.csv` - Phenotype data including cancer subtypes

Both files should be placed in the appropriate data directory as specified in the file paths.

Key Components

1. Data Preprocessing
   - Batch correction normalization
   - Feature selection and dimensionality reduction
   - Train/test splitting

2. Model Implementation
   - SVC with linear kernel
   - XGBoost classifier with multiclass objective
   - Deep neural network with multiple dense layers

3. Evaluation Metrics
   - Accuracy
   - Precision
   - F1 Score 
   - Recall
   - ROC AUC
   - Confusion matrices

Usage

The code is organized in a Jupyter notebook format. To run the analysis:

1. Ensure all dependencies are installed
2. Place data files in the correct directory
3. Run cells sequentially in the notebook

The notebook will output model performance metrics and visualization plots for analysis.

Results

The models achieve varying performance metrics across different approaches:

- MLP with K-best features achieves approximately 82.5% test accuracy
- SVC with PCA-200 features achieves approximately 80% test accuracy
- XGBoost shows competitive performance with around 77.5% accuracy

Authors

Mihika Sonalkar, Sam Reynolds, George Kucera, and Audrey Choi

License

This project is licensed under standard academic use terms for the course MED 264: Principles of Biomedical Informatics offered by the University of California, San Diego.
