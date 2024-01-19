import csv
import numpy as np
import matplotlib.pyplot as plt

names = ["de Casteljau", "Rational de Casteljau", "Volk-Schumaker algorithm",
         "HornBez algorithm", "lader algorihtm", "Wozny-Chudy algorithm",
         "barycentric", "Wang-Ball algorithm"]

data = []
with open("sample_result.csv", 'r') as csv_file:
    csv_reader = csv.reader(csv_file)

    # Read and print each row of the CSV file
    for row in csv_reader:
        data.append(row)
        
data = np.array(data, dtype=float)
data = data.transpose()

plt.figure(figsize=(15,10))
for index, item in enumerate(data):
    plt.plot(np.array([1, 10, 20, 30, 40, 50])*1000, item, label=names[index], color=f'C{index}')
    
plt.legend()
plt.savefig("comparison_with_respect_to_sample.png", bbox_inches='tight')

data = []
with open("degree_result.csv", 'r') as csv_file:
    csv_reader = csv.reader(csv_file)

    # Read and print each row of the CSV file
    for row in csv_reader:
        data.append(row)
data = np.array(data, dtype=float)
data = data.transpose()

plt.figure(figsize=(15,10))
for index, item in enumerate(data):
    plt.plot(np.array([3, 5, 7, 10, 15, 20]), item, label=names[index], color=f'C{index}')
    
plt.legend()
plt.savefig("comparison_with_respect_to_degree.png", bbox_inches='tight')