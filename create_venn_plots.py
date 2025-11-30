#!/usr/bin/env python3
"""
Venn Diagram Generator for Peptide Sequence Predictions
Creates Venn diagrams for comparing predictions from multiple de novo sequencing tools.

Supports:
- 2-3 way Venn diagrams using matplotlib-venn
- 4-6 way Venn diagrams using pyvenn
- Aborts if more than 6 tools are selected

Usage:
    python venn_plot_generator.py <csv_file>
    python venn_plot_generator.py <csv_file> --tools tool1 tool2 tool3
    python venn_plot_generator.py <csv_file> --output custom_output.png
    python venn_plot_generator.py <csv_file> --title "Custom Title"
    python venn_plot_generator.py <csv_file> --dpi 300
    python venn_plot_generator.py <csv_file> --figsize 12 10
    python venn_plot_generator.py <csv_file> --exact  # Use exact string matching instead of mass-based
"""

import argparse
import re
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3, venn2_circles, venn3_circles
import numpy as np
from pathlib import Path
import sys
import warnings
warnings.filterwarnings('ignore')

import calculate_criteria
from constants import aa_dict

# Try to import pyvenn for 4-6 way diagrams
try:
    import venn as pyvenn
    PYVENN_AVAILABLE = True
except ImportError:
    PYVENN_AVAILABLE = False

from adjustText import adjust_text


def read_predictions(csv_file, selected_tools=None):
    """
    Read the CSV file and extract prediction columns.
    
    Parameters
    ----------
    csv_file : str
        Path to the CSV file
    selected_tools : list or None
        List of specific tools to include, or None for all tools
    
    Returns
    -------
    df : pd.DataFrame
        The loaded dataframe
    tool_columns : list
        List of tool names found in the file
    """
    try:
        df = pd.read_csv(csv_file)
    except FileNotFoundError:
        print(f"Error: File '{csv_file}' not found.")
        sys.exit(1)
    except pd.errors.EmptyDataError:
        print(f"Error: File '{csv_file}' is empty.")
        sys.exit(1)
    
    # Find all tool columns (ending with _seq but not groundtruth_seq)
    seq_columns = [col for col in df.columns if col.endswith('_seq')]
    available_tools = [col.replace('_seq', '') for col in seq_columns if col != 'groundtruth_seq']
    
    if not available_tools:
        print("Error: No tool prediction columns found (looking for columns ending with '_seq')")
        sys.exit(1)
    
    if 'groundtruth_seq' not in df.columns:
        print("Error: 'groundtruth_seq' column not found in the CSV file")
        sys.exit(1)
    
    print(f"Available tools in file: {', '.join(available_tools)}")
    
    # Select specific tools if requested
    if selected_tools:
        # Validate that requested tools exist
        missing_tools = [tool for tool in selected_tools if tool not in available_tools]
        if missing_tools:
            print(f"Error: The following tools were not found in the file: {', '.join(missing_tools)}")
            print(f"Available tools are: {', '.join(available_tools)}")
            sys.exit(1)
        tool_columns = selected_tools
        print(f"Using selected tools: {', '.join(tool_columns)}")
    else:
        tool_columns = available_tools
        print(f"Using all {len(tool_columns)} tools")
    
    # Check number of tools
    if len(tool_columns) < 2:
        print("Error: Need at least 2 tools to create a Venn diagram")
        sys.exit(1)
    elif len(tool_columns) > 6:
        print(f"Error: Cannot create Venn diagram for {len(tool_columns)} tools (maximum is 6)")
        print("Please select up to 6 tools using the --tools argument")
        print("Example: --tools tool1 tool2 tool3 tool4 tool5 tool6")
        sys.exit(1)
    
    print(f"Processing {len(df)} peptide sequences")
    
    return df, tool_columns


def get_correct_predictions_sets(df, tool_columns, use_exact=False, aa_dict=None):
    """
    Create sets of indices where each tool made correct predictions.
    Can use either exact string matching or mass-based matching.
    
    Parameters
    ----------
    df : pd.DataFrame
        The dataframe with predictions
    tool_columns : list
        List of tool names
    use_exact : bool
        If True, use exact string matching. If False, use mass-based matching
    aa_dict : Dict[str, float] or None
        Mapping of amino acid tokens to their mass values. If None, uses default.

    Returns
    -------
    sets_dict : dict
        Dictionary mapping tool names to sets of correct prediction indices
    """
    sets_dict = {}

    if use_exact:
        print("Calculating correct predictions (exact string matching)")

    else:
        print("Calculating correct predictions (mass-based matching)")
    
    
    for tool in tool_columns:
        tool_col = f"{tool}_seq"
        # Handle missing values by replacing them with empty strings
        tool_predictions = df[tool_col].fillna('')
        groundtruth = df['groundtruth_seq'].fillna('')
        
        # Find indices where predictions match ground truth
        correct_mask = tool_predictions == groundtruth
        correct_indices = set(df.index[correct_mask].tolist())

        if use_exact:
            # Exact string matching
            for idx, (pred, truth) in enumerate(zip(tool_predictions, groundtruth)):
                if pred and truth and pred == truth:
                    correct_indices.add(idx)

        else:
            # Mass-based matching
            for idx, (pred, truth) in enumerate(zip(tool_predictions, groundtruth)):
                if pred and truth:  # Skip empty sequences
                    # Tokenize peptides
                    pred_tokens = re.split(r'(?<=.)(?=[A-Z])', pred)
                    truth_tokens = re.split(r'(?<=.)(?=[A-Z])', truth)
                    
                    # Use aa_match to determine if they match
                    aa_matches, pep_match = calculate_criteria.aa_match(
                        pred_tokens, truth_tokens, aa_dict
                    )
                    
                    if pep_match:  # Full peptide match based on mass
                        correct_indices.add(idx)
        
        sets_dict[tool] = correct_indices
        
        # Print statistics
        n_correct = len(correct_indices)
        n_total = len(df)
        accuracy = n_correct / n_total * 100 if n_total > 0 else 0
        print(f"  {tool}: {n_correct}/{n_total} correct ({accuracy:.2f}%)")
    
    return sets_dict


def create_venn2(sets_dict, labels, title, output_file, figsize=(8, 6), dpi=150, show_numbers=True):
    """
    Create a 2-way Venn diagram.
    """
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    
    # Get the two sets
    set_a, set_b = [sets_dict[label] for label in labels]
    
    # Create the Venn diagram
    venn = venn2([set_a, set_b], set_labels=labels, ax=ax)
    
    # Add circles for better visibility
    venn2_circles([set_a, set_b], ax=ax, linewidth=1.5)
    
    # Style the diagram
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
    
    # Handle number display
    if show_numbers:
        # Add percentage annotations
        total = len(set_a.union(set_b)) if set_a.union(set_b) else 1
        for subset in ['10', '01', '11']:
            if venn.get_label_by_id(subset):
                count = int(venn.get_label_by_id(subset).get_text())
                percentage = count / total * 100
                venn.get_label_by_id(subset).set_text(f'{count}\n({percentage:.1f}%)')
    else:
        # Hide all region labels
        for subset in ['10', '01', '11']:
            if venn.get_label_by_id(subset):
                venn.get_label_by_id(subset).set_text('')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
    print(f"\nVenn diagram saved to: {output_file}")
    plt.close()


def create_venn3(sets_dict, labels, title, output_file, figsize=(10, 8), dpi=150, show_numbers=True):
    """
    Create a 3-way Venn diagram.
    """
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    
    # Get the three sets
    set_a, set_b, set_c = [sets_dict[label] for label in labels]
    
    # Create the Venn diagram
    venn = venn3([set_a, set_b, set_c], set_labels=labels, ax=ax)
    
    # Add circles for better visibility
    venn3_circles([set_a, set_b, set_c], ax=ax, linewidth=1.5)
    
    # Style the diagram
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
    
    # Handle number display
    if show_numbers:
        # Add percentage annotations
        total = len(set_a.union(set_b).union(set_c)) if set_a.union(set_b).union(set_c) else 1
        for subset in ['100', '010', '001', '110', '101', '011', '111']:
            if venn.get_label_by_id(subset):
                count = int(venn.get_label_by_id(subset).get_text())
                percentage = count / total * 100
                venn.get_label_by_id(subset).set_text(f'{count}\n({percentage:.1f}%)')
    else:
        # Hide all region labels
        for subset in ['100', '010', '001', '110', '101', '011', '111']:
            if venn.get_label_by_id(subset):
                venn.get_label_by_id(subset).set_text('')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
    print(f"\nVenn diagram saved to: {output_file}")
    plt.close()


def create_pyvenn_diagram(sets_dict, labels, title, output_file, figsize=(12, 10), dpi=150, show_numbers=True):
    """
    Create 4, 5, or 6-way Venn diagrams using pyvenn.
    """
    if not PYVENN_AVAILABLE:
        print("\nError: pyvenn is required for 4-6 way Venn diagrams")
        print("Please install it with: pip install pyvenn")
        sys.exit(1)
    
    n_sets = len(labels)
    if n_sets < 4 or n_sets > 6:
        print(f"Error: pyvenn only supports 4-6 sets, got {n_sets}")
        sys.exit(1)
    
    # Convert sets to pyvenn format - it needs a dictionary of labels to sets
    labels_dict = {label: sets_dict[label] for label in labels}
    
    # Set up matplotlib figure with specified size and DPI
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    
    # Create the Venn diagram - use empty format string if hiding numbers
    fmt = "{size}" if show_numbers else ""
    venn_diagram = pyvenn.venn(labels_dict, ax=ax, fmt=fmt, fontsize=8)
    
    # Only apply adjust_text if showing numbers
    if show_numbers:
        # Collect all text objects from the venn diagram for adjustment
        texts = [child for child in ax.get_children() 
                 if isinstance(child, plt.Text) and child.get_text().isdigit()]
        
        # Use adjust_text to prevent overlapping labels
        if texts:
            adjust_text(
                texts,
                ax=ax,
                expand_points=(1.5, 1.5),
                expand_text=(1.2, 1.2),
                force_points=(0.5, 0.5),
                force_text=(0.5, 0.5),
                arrowprops=dict(arrowstyle='-', color='gray', lw=0.5, alpha=0.5)
            )
    
    # Add title
    ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=dpi, bbox_inches='tight')
    print(f"\nVenn diagram saved to: {output_file}")
    plt.close()



def print_intersection_statistics(sets_dict):
    """
    Print basic statistics about intersections.
    """
    tools = list(sets_dict.keys())
    all_correct = set().union(*sets_dict.values())
    
    print("\n" + "="*60)
    print("INTERSECTION STATISTICS")
    print("="*60)
    
    # Overall statistics
    print(f"\nTotal unique peptides correctly predicted by at least one tool: {len(all_correct)}")
    
    # Find peptides predicted by all tools
    if len(tools) > 1:
        all_tools_correct = sets_dict[tools[0]].copy()
        for tool in tools[1:]:
            all_tools_correct = all_tools_correct.intersection(sets_dict[tool])
        print(f"Peptides correctly predicted by ALL tools: {len(all_tools_correct)}")
    
    # Find unique predictions for each tool
    print("\nUnique correct predictions per tool (not in any other tool):")
    for tool in tools:
        unique = sets_dict[tool].copy()
        for other_tool in tools:
            if other_tool != tool:
                unique = unique - sets_dict[other_tool]
        print(f"  {tool}: {len(unique)} unique")
    
    print("="*60)


def main():

    parser = argparse.ArgumentParser(
        description='Generate Venn diagrams from peptide prediction results (max 6 tools)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
    Examples:
    # Use all tools (aborts if more than 6)
    python venn_plot_generator.py my_dataset_benchmark_predictions.csv
    
    # Select specific tools
    python venn_plot_generator.py data.csv --tools tool1 tool2 tool3
    
    # Select up to 6 tools from a file with more tools
    python venn_plot_generator.py data.csv --tools casanovo novor pointnovo pepnet directms specgpt
    
    # Custom output and styling
    python venn_plot_generator.py data.csv --output results/venn.png --dpi 300
    python venn_plot_generator.py data.csv --title "Tool Comparison" --figsize 14 12
            """
    )
    
    parser.add_argument('csv_file', help='Path to the CSV file containing predictions')
    parser.add_argument('--tools', '-t', nargs='+', 
                       help='Specific tools to include (max 6). If not specified, uses all tools.')
    parser.add_argument('--output', '-o', help='Output file path (default: auto-generated based on input)')
    parser.add_argument('--title', help='Custom title for the plot')
    parser.add_argument('--dpi', type=int, default=150, help='DPI for the output image (default: 150)')
    parser.add_argument('--figsize', nargs=2, type=float, default=[12, 10], 
                       help='Figure size as width height (default: 12 10)')
    parser.add_argument('--exact', action='store_true',
                       help='Use exact string matching instead of mass-based matching (default: mass-based)')
    parser.add_argument('--no-numbers', action='store_true',
                       help='Hide the count numbers in the Venn diagram regions')
    
    args = parser.parse_args()
    
    # Validate number of selected tools
    if args.tools and len(args.tools) > 6:
        print(f"Error: Too many tools selected ({len(args.tools)}). Maximum is 6.")
        print(f"You selected: {', '.join(args.tools)}")
        sys.exit(1)
    
    # Determine output file path
    input_path = Path(args.csv_file)
    if args.output:
        output_file = args.output
    else:
        # Auto-generate output filename
        dataset_name = input_path.stem.replace('_benchmark_predictions', '')
        modus = 'string_match' if args.exact else 'mass_match'
        output_file = input_path.parent / f"{dataset_name}_predictions_venn_plot_{modus}.png"



    
    # Determine title (empty if not provided)
    title = args.title if args.title else ""
    
    print(f"\nProcessing: {args.csv_file}")
    print(f"Output will be saved to: {output_file}")
    
    # Read the CSV file
    df, tool_columns = read_predictions(args.csv_file, args.tools)
    
    # Get correct predictions for each tool
    sets_dict = get_correct_predictions_sets(df, tool_columns, use_exact=args.exact, aa_dict=aa_dict)
    
    # Print statistics
    print_intersection_statistics(sets_dict)
    
    # Create Venn diagram based on number of tools
    n_tools = len(tool_columns)
    
    print(f"\nCreating {n_tools}-way Venn diagram...")
    
    if n_tools == 2:
        create_venn2(sets_dict, tool_columns, title, output_file, 
                    figsize=tuple(args.figsize), dpi=args.dpi, show_numbers=not args.no_numbers)
    elif n_tools == 3:
        create_venn3(sets_dict, tool_columns, title, output_file, 
                    figsize=tuple(args.figsize), dpi=args.dpi, show_numbers=not args.no_numbers)
    elif n_tools >= 4 and n_tools <= 6:
        create_pyvenn_diagram(sets_dict, tool_columns, title, output_file,
                            figsize=tuple(args.figsize), dpi=args.dpi, show_numbers=not args.no_numbers)
    
    print("\nDone!")


if __name__ == "__main__":
    main()