import umap
from matplotlib import pyplot as plt
from sklearn.datasets import load_digits


digits = load_digits()

embedding = umap.UMAP().fit_transform(digits.data)

plt.scatter(embedding[:, 0], embedding[:, 1])

a = 1