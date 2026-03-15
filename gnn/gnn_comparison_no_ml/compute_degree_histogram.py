"""
Compute the degree distribution histogram of the dataset (used for the PNA model)
"""
import torch
from torch_geometric.utils import degree

# Load dataset
dataset = torch.load('../processed_boron_dataset_v2.pt')

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
