import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

from scripts import linkageList

def heat_dendrogram(truthJet=None, recluster_jet1=None, recluster_jet2=None, heat_map=True, full_path=False, FigName = None):


	if truthJet:

		# Calculate linkage list and add it to the jet dict
		linkageList.draw_truth(truthJet)

		level_length = [len(entry) for entry in truthJet["tree_ancestors"]]
		max_level = np.max(level_length)

		ancestors1_array = np.asarray([np.concatenate(
			(truthJet["tree_ancestors"][i], -(i + 1) * np.ones((max_level - len(truthJet["tree_ancestors"][i])))))
		                               for i in range(len(truthJet["tree_ancestors"]))])

	else:

		ancestors1 = recluster_jet1["tree_ancestors"]
		level_length = [len(entry) for entry in ancestors1]
		max_level = np.max(level_length)

		ancestors1_array = np.asarray(
			[np.concatenate((ancestors1[i], -(i + 1) * np.ones((max_level - len(ancestors1[i])))))
			 for i in range(len(ancestors1))])

	N_heat = len(ancestors1_array)
	heat_data = np.zeros((N_heat, N_heat))
	neg_entries = np.sum(np.array(ancestors1_array) < 0, axis=1)

	if full_path:
		for i in range(N_heat):
			for j in range(i + 1, N_heat):
				#             print(np.count_nonzero(tree_ancestors_1_array[i]-tree_ancestors_1_array[j]==0))
				heat_data[i, j] = np.absolute(2 * np.count_nonzero((ancestors1_array[i] - ancestors1_array[j]) == 0) -
				                              (2 * max_level - (neg_entries[i] + neg_entries[j])))
				heat_data[j, i] = heat_data[i, j]

	else:
		for i in range(N_heat):
			for j in range(i + 1, N_heat):
				#             print(np.count_nonzero(tree_ancestors_1_array[i]-tree_ancestors_1_array[j]==0))
				heat_data[i, j] = np.absolute(np.count_nonzero((ancestors1_array[i] - ancestors1_array[j]) == 0) -
				                              (max_level - np.minimum(neg_entries[i], neg_entries[j])))
				heat_data[j, i] = heat_data[i, j]

	if truthJet:
		print('---' * 20)
		print('truth heat data ----  alpha row: truth -- alpha column:', 'truth')
		sns.clustermap(
			heat_data,
			row_cluster=True,
			col_cluster=True,
			row_linkage=truthJet["linkage_list"],
			col_linkage=truthJet["linkage_list"],
		)

		if FigName:
			plt.savefig(str(FigName))

		plt.show()

		if recluster_jet1:

			print('---' * 20)
			print('alpha row:', recluster_jet1["algorithm"], '-- alpha column:', 'truth')
			sns.clustermap(
				heat_data,
				row_cluster=True,
				col_cluster=True,
				row_linkage=recluster_jet1["linkage_list"],
				col_linkage=truthJet["linkage_list"],
			)
			plt.show()


	else:
		print('---' * 20)
		print('alpha row:', recluster_jet1["algorithm"], '-- alpha column:', recluster_jet1["algorithm"])
		sns.clustermap(
			heat_data,
			row_cluster=True,
			col_cluster=True,
			row_linkage=recluster_jet1["linkage_list"],
			col_linkage=recluster_jet1["linkage_list"],
		)
		plt.show()

		if recluster_jet2:
			print('---' * 20)
			print('alpha row:', recluster_jet2["algorithm"], '-- alpha column:', recluster_jet1["algorithm"])
			sns.clustermap(heat_data, pivot_kws=None, z_score=None,
			               standard_scale=None, figsize=None, cbar_kws=None, row_cluster=True, col_cluster=True,
			               row_linkage=recluster_jet2["linkage_list"], col_linkage=recluster_jet1["linkage_list"],
			               row_colors=None, col_colors=None, mask=None)
			plt.show()