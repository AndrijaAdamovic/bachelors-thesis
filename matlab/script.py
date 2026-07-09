import csv

filename = "results.csv"

with open(filename, "r") as file:
    reader = csv.reader(file)
    
    # If your CSV has a header row, uncomment the next line:
    # headers = next(reader)

    # Read all rows and convert values to floats
    data = [[float(val) for val in row] for row in reader]

# Transpose rows into columns
columns = list(zip(*data))
num_rows = len(data)

print("--- Column Metrics ---")
for idx, col in enumerate(columns):
    col_avg = sum(col) / num_rows
    col_max = max(col)
    col_min = min(col)
    
    print(f"Column {idx + 1}:")
    print(f"  Average: {col_avg}")
    print(f"  Largest: {col_max}")
    print(f"  Smallest: {col_min}")
    print("-" * 22)
