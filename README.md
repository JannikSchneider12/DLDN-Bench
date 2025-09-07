# Benchmarking & Evaluation of Deep Learning - based De Novo Peptide Sequencing Models

## General process

1. Get the Evaluation Datasets (Benchmark MGF files) and apply your new De Novo Model to get predictions
2. Clone this repo and apply the evaluation to see how your new model performs with respect to current state-of-the-art models on the benchmark datasets
    2.1. Your predictions need to have a specific format:

3. Run the create_result_csv.py which takes your predictions, preprocesses it and aligns it with the provided state-of-the-art models predictions
4. Run the calc_and_plot_precision_coverage.py which calculates the evaluation criteria and generates the precision-coverage plots on the peptide and on the aminoacid level as well as the AUC-score (area under the curve)