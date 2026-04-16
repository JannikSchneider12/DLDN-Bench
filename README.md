# Benchmarking & Evaluation of Deep Learning - based De Novo Peptide Sequencing Models

## General process

1. Get the Evaluation Datasets (Benchmark MGF files) which are uploaded on Zenodo: https://doi.org/10.5281/zenodo.18213583 and apply your new De Novo Model to get predictions

2. Clone this repo to apply the evaluation to see how your new model performs with respect to current state-of-the-art models on the benchmark datasets
    
    - 2.1. Your predictions need to have a specific format:

        - A csv file with at least the following columns:

        | predictions   | scores                       | 
        |---------------|------------------------------|
        | SEQ1          | Score1                       | 
        | SEQ2          | Score2                       |
        | ...           | ...                          | 

        - the predictions column should contain the predicted sequence and the scores should be the models confidence score or assigned probabilities
        - Note on modifications: If your predictions contain modifications that are not supported by our evaluation pipeline or have a different format, please edit the dictionaries in "constants.py" accordingly
        - Note on Indexing: The evaluation uses the positional index, so the order should be the same in your provided prediction file, e.g. entry 7 should contain the prediction for spectrum 7 in the benchmark MGF. If the model that you want to test shuffles the order and uses some identifier instead (e.g. the spectrum title), you need to map this to the correct positional order again such that it aligns with the order of the original benchmark file 

3. Run the create_result_csv.py which takes your predictions, preprocesses it and aligns it with the provided state-of-the-art models predictions

4. Run the calc_and_plot_precision_coverage.py which calculates the evaluation criteria and generates the precision-coverage plots on the peptide and on the aminoacid level as well as the AUC-score (area under the curve)
