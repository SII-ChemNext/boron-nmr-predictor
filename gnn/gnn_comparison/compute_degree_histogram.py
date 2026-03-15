"""
Compute the degree distribution histogram of the dataset (used for the PNA model).
Uses the processed_boron_dataset_v3.pt in the current directory.
"""
import torch
from torch_geometric.utils import degree

# Load dataset
DATA_PATH = 'processed_boron_dataset_v3.pt'
print(f"Loading dataset: {DATA_PATH}")
dataset = torch.load(DATA_PATH, weights_only=False)
print(f"Number of samples: {len(dataset)}")

# Compute maximum degree
max_degree = 0
for data in dataset:
    d = degree(data.edge_index[1], num_nodes=data.num_nodes, dtype=torch.long)
    max_degree = max(max_degree, int(d.max()))

# Compute degree histogram
deg = torch.zeros(max_degree + 1, dtype=torch.long)
for data in dataset:
    d = degree(data.edge_index[1], num_nodes=data.num_nodes, dtype=torch.long)
    deg += torch.bincount(d, minlength=deg.numel())

print(f"Maximum degree: {max_degree}")
print(f"Degree histogram: {deg}")

# Save
torch.save(deg, 'pna_deg_histogram.pt')
print("Saved to pna_deg_histogram.pt")
