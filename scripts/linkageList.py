import sys
import numpy as np
import logging
import pickle
import logging

from scripts.utils import get_logger

logger = get_logger(level=logging.INFO)





def draw_truth(in_jet):
	"""
	Specific function to build the linkage list for the truth jet tree only.
	This is a special case because we build the linkage list from a top down approach.
	(Typically it is built from a bottom-up approach)

	[SciPy linkage list website](https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html)
	Linkage list format: A  (n - 1) by 4 matrix Z is returned. At the i-th iteration, clusters with indices Z[i, 0] and Z[i, 1] are combined to form cluster (n + 1) . A cluster with an index less than n  corresponds to one of the n original observations. The distance between clusters Z[i, 0] and Z[i, 1] is given by Z[i, 2]. The fourth value Z[i, 3] represents the number of original observations in the newly formed cluster.

	Args:
	- input jet dictionary.
	"""

	runTraverse_jet(in_jet, draw_tree=True)

	outers_node_id = in_jet["outers_node_id"]

	N = outers_node_id[-1]
	Nleaves = len(outers_node_id)

	const_list = np.arange(Nleaves)
	temp_outers_node_id = outers_node_id
	idx = np.asarray(outers_node_id)
	temp = []

	sibling_pairs = np.asarray(list(in_jet["parent_child"].values()))[::-1] # All sibling pairs node idxs
	logger.debug(f"sibling_pairs = {sibling_pairs}")

	N_leaves_list = np.ones((Nleaves)) # List that given a node idx, stores for that idx, the number of leaves for the branch below that node.
	linkage_list = [] # linkage list to build heat clustermap visualizations.
	waitlist = {}
	j = 0
	m = 0

	# Start from the last leaf and build the linkage list going backwards
	for i in range(N, 0, -1):
		m += 1
		node_loc = np.where(sibling_pairs == i) # Get location of node in sibling_pairs list
		sibling = sibling_pairs[node_loc[0], 1 - node_loc[1]] # Get sibling node idx

		logger.debug(f"Location of node in sibling_pairs list = {node_loc}")
		logger.debug(f"sibling = {sibling}")

		node_pos = None
		if i in temp_outers_node_id:
			node_pos = np.where(idx == i)[0][0] # This will give the position of i in const_list
			temp_outers_node_id = np.delete(temp_outers_node_id, node_pos) # If i is the idx of a leaf, then delete i from list of leaves

		else:
			const_list = np.append(const_list, Nleaves + j) # Append i as an inner node to the const list
			node_pos = Nleaves + j
			j += 1


		if node_pos != None:

			# If the sibling is a leaf
			if sibling[0] in temp_outers_node_id:

				sibling_pos = np.where(idx == sibling[0]) #Find sibling position in outers list
				N_leaves_list = np.append(N_leaves_list, N_leaves_list[node_pos] + N_leaves_list[sibling_pos[0]])

				linkage_list.append([np.minimum(const_list[node_pos], const_list[sibling_pos[0]][0]),
				                     np.maximum(const_list[node_pos], const_list[sibling_pos[0]][0]), m,
				                     N_leaves_list[-1]])

			# If sibling is an inner node that was stored in waitlist
			elif sibling[0] in waitlist.keys():
				N_leaves_list = np.append(N_leaves_list, N_leaves_list[node_pos] + N_leaves_list[waitlist[sibling[0]]])

				linkage_list.append([np.maximum(const_list[node_pos], const_list[waitlist[sibling[0]]]),
				                     np.minimum(const_list[node_pos], const_list[waitlist[sibling[0]]]), m,
				                     N_leaves_list[-1]])

			# If we do not have the sibling yet, we store the node in waitlist.
			else:
				waitlist[i] = node_pos


		logger.debug(f"linkage_list = {linkage_list}")
		logger.debug(f"waitlist = {waitlist}")


	in_jet["linkage_list"] = np.asarray(linkage_list)








def runTraverse_jet(in_jet, parent=-1, node_id=0, dendrogram=True,
                    ancestors=[], tree_ancestors=[], draw_tree=False):
	'''
	Runs the traverse_jet function.

	Args:
	- in_jet: input jet dictionary.
	- parent: parent node id for the starting point, if -1 then we start from the root of the tree.
	- node_id: node id of the starting node, root node if node_id=0.
	- draw_tree: Bool. Flag that saves "parent_child" and "outers_node_id" lists in the jet dictionary.
	Returns:
	 outers_list, tree_ancestors, ( parent_child, outers_node_id )

	'''

	parent_child_dic = {}
	tree_ancestors = []
	outers_list = []
	outers_node_id = []

	node_id = in_jet['root_id']

	outers_list, tree_ancestors, parent_child, outers_node_id = traverse_jet(
		in_jet,
		parent=parent,
		node_id=node_id,
		outers_list=outers_list,
		dendrogram=dendrogram,
		ancestors=ancestors,
		tree_ancestors=tree_ancestors,
		parent_child_dic=parent_child_dic,
		outers_node_id=outers_node_id,
	)


	in_jet["outers_list"] = outers_list
	in_jet["tree_ancestors"] = tree_ancestors

	if draw_tree:
		in_jet["parent_child"] = parent_child
		in_jet["outers_node_id"] = outers_node_id






def traverse_jet(
		jet,
		parent = None,
		node_id = None,
		outers_list = None,
		dendrogram = True,
		ancestors = None,
		tree_ancestors = [],
		parent_child_dic = None,
		outers_node_id = None,
):
	'''
	Recursive function to traverse the tree and get a list of the leaves
	Args:
		jet: jet dictionary
		node_id: id of the current node
		outers_list: list that stores the momentum of the leaves of the tree
		ancestors: 1D array with all the ancestors from root to leaf for each leaf.
		tree_ancestors: List with one entry for each leaf of the tree, where each entry lists all the ancestor node ids when traversing the tree from the root to the leaf node. (Each entry is one "ancestors" array)
		parent_child_dic: dictionary where each key is a parent node and the values a list of the children
		outers_node_id: list that stores the node id of the outers in the order in which we traverse the tree.

	Returns:
		outers_list, tree_ancestors, parent_child_dic, outers_node_id
	'''

	new_ancestors = None
	if dendrogram:
		new_ancestors = np.copy(ancestors)
		new_ancestors = np.append(new_ancestors, node_id)  # Node idxs in the truth jet dictionary

	# Build outers list
	if jet["tree"][node_id, 0] == -1:
		outers_list.append(jet["content"][node_id])
		outers_node_id.append(node_id)

		if dendrogram:
			tree_ancestors.append(new_ancestors)

	else:
		# Add {parent:[children]} to dic
		parent_child_dic[node_id] = jet["tree"][node_id]

		traverse_jet(
			jet,
			node_id,
			jet["tree"][node_id, 0],
			outers_list,
			ancestors=new_ancestors,
			dendrogram=dendrogram,
			tree_ancestors=tree_ancestors,
			parent_child_dic=parent_child_dic,
			outers_node_id=outers_node_id,
		)

		traverse_jet(
			jet,
			node_id,
			jet["tree"][node_id, 1],
			outers_list,
			ancestors=new_ancestors,
			dendrogram=dendrogram,
			tree_ancestors=tree_ancestors,
			parent_child_dic=parent_child_dic,
			outers_node_id=outers_node_id,
		)

	return outers_list, tree_ancestors, parent_child_dic, outers_node_id