from graphviz import Digraph
import numpy as np
import logging

from scripts.utils import get_logger

logger = get_logger(level=logging.INFO)


def plotBinaryTree(
		jet,
		label = False,
		node_id_in = None,
		pySort = False,
		pTSort = False,
		truthOrder = None
):
	'''
	Plot a tree with nodes and edges using graphviz Digraph
	Args:
	- jet: input jets dictionary
	- label: Bool. If true, then add the label to each node (currently using (py,pz) and Delta as label)
	- node_id_in: index of the node id of each leaf  in the original tree used to get the list of leaves that were
	  subsequently reclustered (currently using the truth level tree order). So, if we sort the node id's of the new algorithm
	  according to node_id_in, we get the list of leaves in the order of the original algorithm.
	- pySort: sort leaves in increasing py.
	- pTSort: sort leaves in increasing pT=abs(py).
	- truthOrder: Bool. If True, then use the truth tree to order the leaves.

	'''

	content = jet["content"]

	arrowsize = "0.1"

	# Dot will be the tree graph
	dot = Digraph(
		graph_attr={"rank": "flow"},
		edge_attr={"arrowsize": arrowsize},
		node_attr={"style": "filled"},
		format="pdf",
	)


	# Create a subgraph to plot all the leaves at the same level.
	# Connect all the leaves with invisible edges to fix the leaves order within different trees.
	leaves = Digraph(edge_attr={"arrowsize": arrowsize, "style": "invis"})

	# Plot all the leaves at the same level
	leaves.attr(rank='same')

	########################
	outers = []
	# Recursive function to build the graph
	def _rec(jet, parent, node_id):

		# Add a label to each node
		if label:
			if "deltas" in jet.keys() and jet["deltas"][node_id] >= 0.:
				node_label = "py:%0.1f\n pz:%0.1f\n &#916;:%0.2f " % (jet["content"][node_id][0],
				                                                      jet["content"][node_id][1],
				                                                      jet["deltas"][node_id])
			else:
				node_label = "py:%0.1f\n pz:%0.1f" % (jet["content"][node_id][0],
				                                      jet["content"][node_id][1])

		else:
			node_label = ''""


		# Define the subgraph for each recursive call
		sub = Digraph(
			node_attr={"fixedsize": "true",
			           "label": str(node_label),
			           "height": "0.1",
			           "width": "0.1",
			           "style": "filled",
			           "fontsize": "22.0", },
			edge_attr={"arrowsize": arrowsize,
			           "minlen": "0.1"},
		)

		# Nodes size and color
		size = "1.25"
		node_color = "lightblue" if jet["tree"][node_id, 0] == -1 else "wheat"

		# Add node
		sub.node("%d" % node_id,
		         width=size,
		         height=size,
		         shape="circle",
		         color=node_color)


		# Add leaves to outers list
		if jet["tree"][node_id, 0] == -1:
			outers.append(node_id)

		# Add subgraph to main graph
		dot.subgraph(sub)


		# Draw from root to leaves. Connect child to parent node. (1st entry is parent, 2nd entry is child)
		if parent >= 0:
			dot.edge("%d" % parent,
			         "%d" % node_id,
			         color="black")

		# Recursive calls
		if jet["tree"][node_id, 0] != -1:

			_rec(jet,
			     node_id,
			     jet["tree"][node_id, 0],
			     )

			_rec(jet,
			     node_id,
			     jet["tree"][node_id, 1],
			     )


	########################
	# Run the recursive function
	_rec(jet, -1, jet["root_id"])


	# Sort the leaves to match the order in which they are accessed when traverseing a tree from some other clustering algorithm (or truth jet).
	# The order is in node_id_in.
	if node_id_in:
		if not truthOrder:
			outers = [outers[k] for k in node_id_in]
		else:
			new_idx_list = list(zip(outers, node_id_in))
			logger.debug(f"new_idx_list before sorting=  {new_idx_list}")
			new_idx_list = sorted(new_idx_list,
			                      key=lambda x: x[1])  # Sort according to node_id_in.
			logger.debug(f"new_idx_list after sorting=  {new_idx_list}")
			outers = [x for (x, y) in new_idx_list]  # List the node ids in the new order.


	# Sort the leaves in increasing py. (Also commnent/uncomment ptList line below to use absolute value)
	if pySort:
		ptList = [jet["content"][node_id][0] for node_id in outers]
	elif pTSort:
		ptList=[np.absolute(jet["content"][node_id][0]) for node_id in outers]

	if pySort or pTSort:
		new_idx_list = list(zip(outers, ptList))
		logger.debug(f"new_idx_list before sorting=  {new_idx_list}")
		new_idx_list = sorted(new_idx_list,
		                      key=lambda x: x[1])
		logger.debug(f"new_idx_list after sorting=  {new_idx_list}")

		outers = [x for (x, y) in new_idx_list]  # List the node ids in the new order.
		logger.debug(f"outers after= {outers}")


	# Add leaf nodes to the leaves digraph in the new order
	size = "0.8"
	node_color = "lightblue"

	sub_leaf = Digraph(
		node_attr={"fixedsize": "true",
		           "height": "0.1",
		           "width": "0.1",
		           "style": "filled"},
		edge_attr={"arrowsize": arrowsize},
	)

	# Add the 1st leaf
	sub_leaf.node("%d" % outers[0],
	              width=size,
	              height=size,
	              shape="circle",
	              color=node_color,
	              fontsize="17.0")

	# Add the remaining leaves
	for j in range(len(outers) - 1):
		sub_leaf.node("%d" % outers[j + 1],
		              width=size,
		              height=size,
		              shape="circle",
		              color=node_color,
		              fontsize="17.0")

		# Add subgraph to leaves graph
		leaves.subgraph(sub_leaf)

		# Connect all the leaves  sequentially (with invisible edges) so that they are all plotted at the same level
		leaves.edge("%d" % outers[j],
		            "%d" % outers[j + 1],
		            color="black")

	# Add leaves to main graph
	dot.subgraph(leaves)


	return dot


def visualizeTreePair(
		in_jet1,
		in_jet2,
		pySort = False,
		pTSort = False,
		truthOrder=True,
		label=True
):
	'''
	Call plotBinaryTree function to create a representation of the jet tree with graphviz Digraph.
	Compares 2 trees (top/bottom). It loads the 2 jet dictionaries.

	Args:
	- in_jet1: jet 1.
	- in_jet2: jet 2.
	- truthOrder: if True, then the leaves are ordered according to the truth jet. If False, then they are ordered according to the other jet.
	- pySort: sort leaves in increasing py.
	- pTSort: sort leaves in increasing pT=abs(py).
	- label: if True, then add labels with info to each node.

	Note:
	- node_id: index of the node id of each leaf  in the clustering algorithm used to get the list of leaves that were
	  subsequently reclustered (currently using the truth level tree order). So, if we sort the node id's of the new algorithm
	  according to node_id_in, we get the list of leaves in the order of the original algorithm.

	'''

	if truthOrder:
		if in_jet1["algorithm"] == "truth":
			jetTop = in_jet2
			jetBottom = in_jet1
		else:
			jetTop = in_jet1
			jetBottom = in_jet2

		node_id = jetTop["node_id"]

	else:
		if in_jet1["algorithm"] == "truth":
			jetTop = in_jet1
			jetBottom = in_jet2
		else:
			jetTop = in_jet2
			jetBottom = in_jet1

		node_id = jetBottom["node_id"]

	tree1 = plotBinaryTree(
		jetTop,
		label = label,
		node_id_in = node_id,
		pySort = pySort,
		pTSort = pTSort,
		truthOrder = truthOrder,
	)

	tree2 = plotBinaryTree(
		jetBottom,
		label = label,
		pySort = pySort,
		pTSort = pTSort,
	)

	return tree1, tree2