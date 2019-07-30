import sys

import numpy as np
import logging
import pickle


def draw_truth(in_jet):
	"""
	Build linkage list for the truth jet tree.

	Args:
	- input jet dictionary with the tree, content arrays, etc.


	"""

	runTraverse_jet(in_jet, draw_tree=True)

	outers_node_id = in_jet["outers_node_id"]

	N = outers_node_id[-1]
	Nleaves = len(outers_node_id)
	const_list = np.arange(Nleaves)
	temp_outers_node_id = outers_node_id
	idx = np.asarray(outers_node_id)
	temp = []
	sibling_pairs = np.asarray(list(in_jet["parent_child"].values()))[::-1]
	N_leaves_list = np.ones((Nleaves))
	linkage_list = []
	waitlist = {}
	j = 0
	m = 0
	for i in range(N, 0, -1):
		m += 1
		node = np.where(sibling_pairs == i)
		sibling = sibling_pairs[node[0], 1 - node[1]]
		#         print('node =', node)
		#         print('sibling = ', sibling)
		node_pos = None

		if i in temp_outers_node_id:
			node_pos = np.where(idx == i)[0][0]
			temp_outers_node_id = np.delete(temp_outers_node_id, node_pos)

		else:
			const_list = np.append(const_list, Nleaves + j)
			node_pos = Nleaves + j
			j += 1

		if node_pos != None:
			#             print('node_pos=',node_pos)

			if sibling[0] in temp_outers_node_id:
				sibling_pos = np.where(idx == sibling[0])
				N_leaves_list = np.append(N_leaves_list, N_leaves_list[node_pos] + N_leaves_list[sibling_pos[0]])

				linkage_list.append([np.maximum(const_list[node_pos], const_list[sibling_pos[0]][0]),
				                     np.minimum(const_list[node_pos], const_list[sibling_pos[0]][0]), m,
				                     N_leaves_list[-1]])


			elif sibling[0] in waitlist.keys():
				N_leaves_list = np.append(N_leaves_list, N_leaves_list[node_pos] + N_leaves_list[waitlist[sibling[0]]])

				linkage_list.append([np.maximum(const_list[node_pos], const_list[waitlist[sibling[0]]]),
				                     np.minimum(const_list[node_pos], const_list[waitlist[sibling[0]]]), m,
				                     N_leaves_list[-1]])


			else:
				waitlist[i] = node_pos

	#             print('linkage_list = ', linkage_list)

	#     print('waitlist =',waitlist)

	in_jet["linkage_list"] = np.asarray(linkage_list)


#     return np.asarray(linkage_list)




def runTraverse_jet(in_jet, parent=-1, node_id=0, dendrogram=True,
                    ancestors=[], tree_ancestors=[], draw_tree=False):
	'''
	Run traverse_jet function.

	Args:
	- parent: parent node id for the starting point, if -1 then we start from the root of the tree.
	- node_id: node id of the starting node, root node if node_id=0.
	- draw_tree: Flag the returns extra info when reconstructing the truth level tree.

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

	#     print('len(outers_list) =', len(outers_list))

	in_jet["outers_list"] = outers_list
	in_jet["tree_ancestors"] = tree_ancestors

	if draw_tree:
		in_jet["parent_child"] = parent_child
		in_jet["outers_node_id"] = outers_node_id
#         return outers_list, tree_ancestors, parent_child, outers_node_id

#     else:
#         return outers_list, tree_ancestors





def traverse_jet(
		jet, parent=None,
		node_id=None,
		outers_list=None,
		dendrogram=True,
		ancestors=None,
		tree_ancestors=[],
		parent_child_dic=None,
		outers_node_id=None,
):
	'''
	Recursive function to traverse the tree and get a list of the leaves
	Args:
		jet: jet dictionary
		node_id: id of the current node
		outers_list: list that stores the 4-momentum of the leaves
		ancestors: 1D array with all the ancestors from root to leaf for each  a leaf.
		tree_ancestors: 2D array with all the ancestors from root to leaf where each row corresponds to a leaf.
		parent_child_dic: dictionary where each key is a parent node and the values a list of the children
		outers_node_id: list that stores the node id of the outers in the order in which we traverse the tree.

	Returns:
		outers_list, tree_ancestors, parent_child_dic, outers_node_id
	'''

	new_ancestors = None
	if dendrogram:
		new_ancestors = np.copy(ancestors)
		#         print(' ancestors before=', ancestors)
		#         new_ancestors = np.append(new_ancestors,id) # Node ids in terms of the reclustered jet dictionary
		new_ancestors = np.append(new_ancestors, node_id)  # Node ids in terms of the truth jet dictionary

	if jet["tree"][node_id, 0] == -1:
		outers_list.append(jet["content"][node_id])
		outers_node_id.append(node_id)

		if dendrogram:
			tree_ancestors.append(new_ancestors)
	#             print('tree_ancestors=', tree_ancestors)

	else:
		parent_child_dic[node_id] = jet["tree"][node_id]
		#         parent_child_dic.append([node_id,jet["tree"][node_id]])

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