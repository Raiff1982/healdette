import json
import matplotlib.pyplot as plt
import numpy as np

# Read the JSON file
with open('output/codette_antibody_designs_20250912_150658.json', 'r') as f:
    data = json.load(f)

# Extract sequence lengths and scores
lengths = []
scores = []
categories = []

for binder in data['personalized_binders']:
    length = len(binder['sequence'])
    lengths.append(length)
    scores.append(binder['personalization_score'])
    categories.append('Short' if length < 100 else 'Medium' if length < 200 else 'Long')

# Create a figure with two subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

# Plot 1: Sequence Length vs Score
colors = ['red' if c == 'Short' else 'blue' if c == 'Medium' else 'green' for c in categories]
ax1.scatter(lengths, scores, c=colors, s=100)
ax1.set_xlabel('Sequence Length (amino acids)')
ax1.set_ylabel('Personalization Score')
ax1.set_title('Sequence Length vs. Personalization Score')

# Add legend
from matplotlib.lines import Line2D
legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor=c, label=l, markersize=10)
                  for c, l in zip(['red', 'blue', 'green'], ['Short', 'Medium', 'Long'])]
ax1.legend(handles=legend_elements)

# Plot 2: Score Distribution
ax2.bar(['Short', 'Medium', 'Long'], 
        [np.mean([s for s, c in zip(scores, categories) if c == cat]) 
         for cat in ['Short', 'Medium', 'Long']])
ax2.set_ylabel('Average Personalization Score')
ax2.set_title('Average Score by Sequence Category')

# Adjust layout and save
plt.tight_layout()
plt.savefig('output/sequence_analysis.png', dpi=300, bbox_inches='tight')
plt.close()

# Print summary statistics
print("\nSequence Analysis Summary:")
print("-" * 50)
print(f"Total sequences: {len(lengths)}")
print("\nBy Category:")
for category in ['Short', 'Medium', 'Long']:
    cat_scores = [s for s, c in zip(scores, categories) if c == category]
    cat_lengths = [l for l, c in zip(lengths, categories) if c == category]
    if cat_scores:
        print(f"\n{category} sequences:")
        print(f"Count: {len(cat_scores)}")
        print(f"Length range: {min(cat_lengths)}-{max(cat_lengths)} amino acids")
        print(f"Average score: {np.mean(cat_scores):.3f}")
        print(f"Score range: {min(cat_scores):.3f}-{max(cat_scores):.3f}")