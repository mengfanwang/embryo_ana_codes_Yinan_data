import umap
from matplotlib import pyplot as plt
from sklearn.datasets import load_digits
from scipy.io import loadmat, savemat
import numpy as np

 # amat format:  'node id', 'type', 'x', 'y', 'z', 'radius', 'parent node id', 'timepoint', 'tag', 'lineage id'
data = loadmat("/home/mengfan/ForExecute/embryo_ana_codes_Yinan_data/investigation/tracks_stitch_v4_1000.mat")
tracks = data['tracks']
tracks[:, 0] -= 1
tracks[:, 6] -= 1

tracks_last = tracks[tracks[:, 7] == 1000, :]
cnt = 0
locations_list = {}
for ii in range(tracks_last.shape[0]):
    locations = np.empty((1001, 3))
    locations[:] = np.nan
    node = int(tracks_last[ii, 0])
    lineage = tracks[node, 9]
    print('lineage: ', lineage)
    while node >= 0:
        tt = int(tracks[node, 7])
        locations[tt, :] = tracks[node, 2:5]
        node = int(tracks[node, 6])
        # print(node)
    if tt != 0:
        print('lineage not start from zero.')
        assert lineage > 8189
    else:
        locations_list[cnt] = locations
        cnt += 1

lineage_num = len(locations_list)
# for ii in range(lineage_num):
#     for tt in range(1001):
#         if np.isnan(locations_list[ii][tt, 0]):
#             print(ii, tt)

for ii in range(lineage_num):
    for tt in range(1001):
        if np.isnan(locations_list[ii][tt, 0]):
            assert 0 < tt < 1000
            if not np.isnan(locations_list[ii][tt-1, 0]) and not np.isnan(locations_list[ii][tt+1, 0]):
                locations_list[ii][tt, :] = (locations_list[ii][tt-1, :] + locations_list[ii][tt+1, :]) / 2
            elif np.isnan(locations_list[ii][tt-1, 0]):
                assert tt-1 != 0
                locations_list[ii][tt-1, :] = locations_list[ii][tt-2, :] * 2/3 + locations_list[ii][tt+1, :] / 3
                locations_list[ii][tt, :] = locations_list[ii][tt-2, :] / 3 + locations_list[ii][tt+1, :] * 2 / 3
            else:
                assert tt+1 != 1000
                locations_list[ii][tt, :] = locations_list[ii][tt-1, :] * 2/3 + locations_list[ii][tt+2, :] / 3
                locations_list[ii][tt+1, :] = locations_list[ii][tt-1, :] / 3 + locations_list[ii][tt+2, :] * 2 / 3

umap_input = np.zeros((lineage_num, 3003))
for ii in range(lineage_num):
    umap_input[ii, :] = locations_list[ii].reshape((1, 3003))

embedding = umap.UMAP().fit_transform(umap_input)

plt.scatter(embedding[:, 0], embedding[:, 1])
plt.show()
a = 1