import os
import csv
import matplotlib.pyplot as plt

def extract_vina_scores(filename):
    scores = []
    with open(filename, 'r') as f:
        for line in f:
            if "REMARK VINA RESULT:" in line:
                score = float(line.split()[3])
                scores.append(score)
                break  # Exit the loop after the first score is found
    return scores

smiles_file = 'naringenin_analogues.txt'
smiles_list = []

with open(smiles_file, 'r') as file:
    reader = csv.reader(file)
    entries = next(reader)  # Read the first row containing the entries

    for i, entry in enumerate(entries):
        index = str(i)  # Convert the integer index to string format
        smiles_list.append(entry.strip().strip("'"))  # Remove leading/trailing spaces and quotation marks if needed

output_folder = "./output/"
file_prefix = "test"

results = []

mevalonate_score = 0

for filename in os.listdir(output_folder):
    if filename.endswith(".txt") and filename.startswith(file_prefix):
        full_path = os.path.join(output_folder, filename)
        scores = extract_vina_scores(full_path)
        if scores:
            avg_score = scores[0]
            test_index = int(filename.lstrip(file_prefix).rstrip(".txt"))
            results.append((test_index, avg_score))

# Plotting
x = [result[0] for result in results]
y = [result[1] for result in results]
print(results)

# Plotting
plt.figure()
plt.plot(x, y, 'o', markeredgewidth=2)
plt.scatter([-200], [-6.8], color='red', marker='o', label='Control 1')
plt.scatter([-100], [-6.8], color='k', marker='o', label='Control 2')
plt.xlabel("Test Index")
plt.ylabel("Affinity Score (kcal/mol)")
plt.title("Docking Score of Naringenin Analogues vs. Control Ligand to PBP2a Active Site")
plt.grid()
plt.legend()

plt.show()

# Creating CSV file
csv_file = "results.csv"
header = ["Index", "SMILES String", "Best Score"]
rows = []

# Sort the results list by index
sorted_results = sorted(results, key=lambda x: x[0])

# Create rows for CSV
rows = []
for result in sorted_results:
    test_index = result[0]
    smiles_string = smiles_list[test_index]  # Adjust index to match list position
    best_score = result[1]
    rows.append([test_index, smiles_string, best_score])

# Write sorted rows to the CSV file
with open(csv_file, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(header)
    writer.writerows(rows)

print(f"CSV file '{csv_file}' created!")
