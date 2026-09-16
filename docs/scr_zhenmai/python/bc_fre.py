from collections import Counter
import pandas as pd
import argparse

ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input", required=True, help="input fq file")
ap.add_argument("-o", "--output", required=True, help="output txt file")

args = vars(ap.parse_args())

input_file = args["input"]
output_file = args["output"]

def analyze_file(filename):
    try:
        with open(filename, 'r') as f:
            # Read lines and strip whitespace/newlines
            lines = [line.strip() for line in f if line.strip()]
        
        # Calculate frequencies
        counter = Counter(lines)
        
        # Create a DataFrame for structured data
        df = pd.DataFrame(counter.items(), columns=['Record', 'Count'])
        
        # Calculate relative frequency (percentage)
        total_records = len(lines)
        df['Frequency (%)'] = (df['Count'] / total_records) * 100
        
        # Sort by count in descending order
        df = df.sort_values(by='Count', ascending=False).reset_index(drop=True)
        
        return df, total_records

    except FileNotFoundError:
        print(f"File {filename} not found.")
        return None, 0

# Usage
df_results, total_count = analyze_file(input_file)

if df_results is not None:
    print(f"Total Records: {total_count}")
    print(f"Unique Records: {len(df_results)}")
    print("\nTop 10 Most Frequent Records:")
    print(df_results.head(10))
    
    # Save to CSV if needed
    df_results.to_csv(output_file, index=False)