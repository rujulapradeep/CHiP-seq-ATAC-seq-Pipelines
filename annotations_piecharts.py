def install_package(package, manager='pip'):
    try:
        import importlib
        importlib.import_module(package)
    except ImportError:
        print(f"{package} is not installed. Installing...")
        if manager == 'pip':
            subprocess.check_call([sys.executable, "-m", "pip", "install", package])
        elif manager == 'conda':
            subprocess.check_call(["conda", "install", "-c", "conda-forge", "-y", package])

def check_dependencies():
    required_packages = ['pandas', 'matplotlib', 'os', 'sys', 'subprocess', 'shutil']
    for package in required_packages:
        install_package(package)
    
    if not shutil.which('pip'):
        print("pip is not installed. Please install pip.")
        sys.exit(1)
    if not shutil.which('conda'):
        print("conda is not installed. Please install conda.")
        sys.exit(1)

#function to sort annotations
def categorize_annotation(annotation):
    if 'Promoter' in annotation:
        return 'Promoters'
    elif 'Intron' in annotation or 'Exon' in annotation:
        return 'Gene Bodies'
    elif 'Distal Intergenic' in annotation:
        return 'Intergenic Regions'
    else:
        return 'Others'

    #function to create piechart
def create_piechart(input_file, group, output_dir, is_female):
    # Read the BED file into a DataFrame
    df = pd.read_csv(input_file, sep='\t', header=None, names=[
        'seqnames', 'chr', 'start', 'end', 'width', 'strand', 'name', 'score',
        'thick.start', 'thick.end', 'thick.width', 'annotation', 'geneChr',
        'geneStart', 'geneEnd', 'geneLength', 'geneStrand', 'geneId',
        'distanceToTSS', 'ENSEMBL', 'SYMBOL', 'GENENAME'
    ])

    # Apply the categorization function to the 'annotation' column
    df['Category'] = df['annotation'].apply(categorize_annotation)

    # Count the occurrences of each category
    category_counts = df['Category'].value_counts()

    # Define colors for each category
    category_colors = {
        'Promoters': '#FF69B4',
        'Gene Bodies': '#800080',
        'Intergenic Regions': '#4169E1',
        'Others': '#00FF7F'
    }

    # Create a DataFrame with only the selected categories
    selected_categories = ['Promoters', 'Gene Bodies', 'Intergenic Regions']
    selected_df = df[df['Category'].isin(selected_categories)]

    # Count the occurrences of selected categories
    selected_category_counts = selected_df['Category'].value_counts()

    # Plot a pie chart with specified colors
    plt.figure(figsize=(8, 8))
    plt.pie(selected_category_counts, labels=selected_category_counts.index, autopct='%1.1f%%', startangle=140, colors=[category_colors[cat] for cat in selected_category_counts.index])
    plt.axis('equal')

    # Determine whether the pie chart is for female or male
    gender = "female" if is_female else "male"

    # Create folder for the group if it doesn't exist
    group_dir = os.path.join(output_dir, group + "_annotations_piechart")
    os.makedirs(group_dir, exist_ok=True)

    # Save the pie chart to a file
    output_filename = os.path.join(group_dir, f"{group}_{gender}_piechart.png")
    plt.savefig(output_filename)
    plt.close()  # Close the plot to release resources

    # Print checkpoint
    print(f"Generated pie chart for {gender} group {group} in {output_filename}")


def run(input_file, output_dir):
    # Read input file to get data and group information
    df = pd.read_csv(input_file, sep='\t', header=None, names=['female_pathway', 'male_pathway', 'group'])

    # Process each row in the input file
    for index, row in df.iterrows():
        create_piechart(row['female_pathway'], row['group'], output_dir, is_female=True)
        create_piechart(row['male_pathway'], row['group'], output_dir, is_female=False)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python annotations_piecharts.py input_file output_directory")
        sys.exit(1)

    input_file = sys.argv[1]
    output_dir = sys.argv[2]

    # Check and install dependencies if needed
    check_dependencies()
    import os
    import sys
    import subprocess
    import pandas as pd
    import matplotlib.pyplot as plt
    import shutil
    # Process input file
    run(input_file, output_dir)