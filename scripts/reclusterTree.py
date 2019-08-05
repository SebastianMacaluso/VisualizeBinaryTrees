import sys
import numpy as np
import logging
import pickle
import itertools

from scripts.utils import get_logger

logger = get_logger(level=logging.INFO)



def recluster(input_jet, alpha=None, save=True):
  """
  Uses helper functions to get the leaves of an  input jet, recluster them following some algorithm determined by the value of alpha,
   create the new tree for the chosen algorithm, make a jet dictionary and save it.

   Create a dictionary with all the jet tree info
  - jet["root_id"]: root node id of the tree
  - jet["content"]: list with the tree nodes (particles) momentum vectors. For the ToyJetsShower we consider a 2D model,
    so we have (py,pz), with pz the direction of the beam axis
  - jet["tree"]: list with the tree structure. Each entry contains a list with the [left,right] children of a node.
    If [-1,-1] then the node is a leaf.

  New features added to the tree:
  - jet["tree_ancestors"]: List with one entry for each leaf of the tree, where each entry lists all the ancestor node ids
    when traversing the tree from the root to the leaf node.
  - jet["linkage_list"]: linkage list to build heat clustermap visualizations.
  - jet["Nconst"]: Number of leaves of the tree.
  - jet["algorithm"]: Algorithm to generate the tree structure, e.g. truth, kt, antikt, CA.

  Args:
  - input_jet: any jet dictionary with the clustering history.
  - alpha: defines the clustering algorithm. alpha={-1,0,1} defines the {anti-kt, CA and kt} algorithms respectively.
  - save: if true, save the reclustered jet dictionary

  Returns:
    jet dictionary
  """


  def _rec(jet, parent, node_id, outers_list):
    """
    Recursive function to get a list of the tree leaves
    """
    if jet["tree"][node_id, 0] == -1:

      outers_list.append(jet["content"][node_id])

    else:
      _rec(
        jet,
        node_id,
        jet["tree"][node_id, 0],
        outers_list,
      )

      _rec(
        jet,
        node_id,
        jet["tree"][node_id, 1],
        outers_list,
      )

    return outers_list


  outers = []

  # Get constituents list (leaves)
  jet_const = np.asarray(
    _rec(
    input_jet,
    -1,
    input_jet["root_id"],
    outers,
  )
  )

  # Run the kt, CA or antikt clustering algorithms
  raw_tree, \
  idx, \
  jet_content, \
  root_node, \
  Nconst, \
  N_leaves_list, \
  linkage_list = ktAntiktCA(jet_const, alpha=alpha)


  # Build the reclustered tree
  tree, \
  content, \
  node_id, \
  tree_ancestors = _traverse(root_node,
                             jet_content,
                             tree_dic=raw_tree,
                             Nleaves=Nconst,
                             )


  # Create jet dictionary with tree features
  jet = {}
  jet["root_id"] = 0
  jet["tree"] = np.asarray(tree).reshape(-1, 2)
  jet["content"] = np.asarray([np.asarray(c) for c in content]).reshape(-1, 2)
  jet["linkage_list"]=linkage_list
  jet["node_id"]=node_id
  jet["tree_ancestors"]=tree_ancestors
  jet["Nconst"]=Nconst
  jet["algorithm"]=alpha


  # Save reclustered tree
  if save:
    out_dir = "data/"
    # print("input_jet[name]=",input_jet["name"])

    algo = str(input_jet["name"]) + '_' + str(alpha)
    out_filename = out_dir + str(algo) + '.pkl'
    logger.info(f"Output jet filename = {out_filename}")
    with open(out_filename, "wb") as f:
      pickle.dump(jet, f, protocol=2)


  return jet






def ktAntiktCA(const_list, alpha=None):
  """
  Runs the dijMinPair function level by level starting from the list of constituents (leaves) until we reach the root of the tree.
  Note: - We refer to both leaves and inner nodes as pseudojets.

  Args:
      - const_list: jet constituents (i.e. the leaves of the tree)
      - alpha: defines the clustering algorithm. alpha={-1,0,1} defines the {anti-kt, CA and kt} algorithms respectively.

  Returns:
      Note:
         Const_list: nodes list after deleting the constituents that are merged and adding the new pseudojet in each level.
          So this should only have the root of the tree at the end.

      - tree_dic: dictionary that has the node id of a parent as a key and a list with the id of the 2 children as the values
      - idx: array that stores the node id
       (the node id determines the location of the momentum vector of a pseudojet in the jet_content array)
        of the pseudojets that are in the current const_list array. It has the same elements as the const_list (they get updated
        level by level).
      - jet_content: array with the momentum of all the nodes of the jet tree (both leaves and inners).
      - root_node: root node id
      - Nconst: Number of leaves of the jet
      - N_leaves_list: List that given a node idx, stores for that idx, the number of leaves for the branch below that node.
      - linkage_list: linkage list to build heat clustermap visualizations.
  """

  Nconst = len(const_list)

  root_node = 2 * Nconst - 2
  logger.debug(f"Root node = (N constituents + N parent) = {root_node}")

  idx = np.arange(Nconst)

  # List that given a node idx, stores for that idx, the number of leaves for the branch below that node.
  # It is initialized only with the tree leaves
  N_leaves_list = np.ones((Nconst))


  linkage_list = []

  dij_hist = []
  tree_dic = {}
  const_list = np.asarray(const_list)
  jet_content = const_list

  for j in range(len(const_list) - 1):
    const_list, \
    dij_hist, \
    tree_dic, \
    idx, \
    jet_content, \
    N_leaves_list, \
    linkage_list = dijMinPair(
      const_list,
      dij_hist,
      tree_dic,
      jet_content,
      idx,
      alpha=alpha,
      Nconst=Nconst,
      Nparent=j,
      N_leaves_list=N_leaves_list,
      linkage_list=linkage_list,
    )

  return tree_dic, idx, jet_content, root_node, Nconst, N_leaves_list, linkage_list







def dijMinPair(
    const_list,
    var_dij_history,
    tree_dic,
    jet_content,
    idx,
    alpha=None,
    Nconst=None,
    Nparent=None,
    N_leaves_list=None,
    linkage_list=None,
):
    """
    -Calculate all d_ij distance (from the generalized kt jet clustering algorithms) between all possible pair of constituents at a certain level and get the minimum.
    -Update the constituents list by deleting the constituents that are merged and adding the new pseudojet
    (We refer to both leaves and inner nodes as pseudojets.)

    Args:
        - const_list: array with the constituents momentum list for the current level (i.e. deleting the constituents that are merged and
          adding the new pseudojet from merging them)
        - var_dij_history: list with all the previous min{d_ij}
        - tree_dic: dictionary that has the node id of a parent as a key and a list with the id of the 2 children as the values
        - jet_content: array with the momentum of all the nodes of the jet tree (both leaves and inners) after adding one
          more level in the clustering.
          We add a new node each time we cluster 2 pseudojets
        - idx: array that stores the node id (the node id determines the location of the momentum of a pseudojet in the jet_content array)
          of the pseudojets that are in the current const_list array. It has the same number of elements as the const_list (they get updated
          level by level).
        - alpha: defines the clustering algorithm. alpha={-1,0,1} defines the {anti-kt, CA and kt} algorithms respectively.
        - Nconst: Number of leaves
        - Nparent: index of each parent added to the tree_dic.
        - N_leaves_list
        - linkage_list: linkage list to build heat clustermap visualizations.
          [SciPy linkage list website](https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html)
          Linkage list format: A  (n - 1) by 4 matrix Z is returned. At the i-th iteration, clusters with indices Z[i, 0] and Z[i, 1] are combined to form cluster (n + 1) . A cluster with an index less than n  corresponds to one of the n original observations. The distance between clusters Z[i, 0] and Z[i, 1] is given by Z[i, 2]. The fourth value Z[i, 3] represents the number of original observations in the newly formed cluster.

    Returns:
        - new_list: new const_list after deleting the constituents that are merged and adding the new pseudojet in the current level.
        - var_dij_history
        - tree_dic
        - idx
        - jet_content
        -  N_leaves_list
        - linkage_list

    """

    # Get all possible pairings
    pairs = np.asarray(list(itertools.combinations(np.arange(len(const_list)), 2)))

    const_list_pt = np.absolute([element[0] for element in const_list])
    logger.debug(f"const_list_pt = {const_list_pt}")


    # Get all dij at current level: dij=min(pTi^(2\alpha),pTj^(2\alpha)) * [arccos((pi.pj)/|pi|*|pj|)]^2
    epsilon=1e-6 #For numerical stability
    dij_list = [(np.sort((const_list_pt[pairs][k]) ** (2 * alpha))[0] * \
                 (np.arccos(np.dot(const_list[pairs][k][0],const_list[pairs][k][1])/
                            (epsilon + np.linalg.norm(const_list[pairs][k][0]) * np.linalg.norm(const_list[pairs][k][1]))\
                            )) ** 2, k)\
                for k in range(len(const_list[pairs]))]

    logger.debug(f"dij_list = {dij_list}")

    cos_arg=(np.sum([np.count_nonzero(np.absolute(np.sum(const_list[pairs][k][0] * const_list[pairs][k][1]) /
                            (np.sqrt(np.sum(const_list[pairs][k][0] ** 2)) * np.sqrt(
                              np.sum(const_list[pairs][k][1] ** 2))))> 1) for k in range(len(const_list[pairs]))]))
    logger.debug(f"Cos arg > 1? = {cos_arg}")

    cosines=[np.absolute(np.sum(const_list[pairs][k][0] * const_list[pairs][k][1]) /
                            (np.sqrt(np.sum(const_list[pairs][k][0] ** 2)) * np.sqrt(
                              np.sum(const_list[pairs][k][1] ** 2)))) for k in range(len(const_list[pairs]))]
    logger.debug(f"pos,value = {[(i,cosines[i]) for i in range(len(cosines)) if np.absolute(cosines[i])<0.99]}")



    # Get pair index (in pairs list) with min dij
    min_tuple = sorted(dij_list, key=lambda x: x[0])[0]
    min_pair = min_tuple[1]
    logger.debug(f"min_pair= {pairs[min_pair]}")


    # List that given a node idx, stores for that idx, the number of leaves for the branch below that node.
    N_leaves_list = np.concatenate(
      (N_leaves_list, [N_leaves_list[idx[pairs[min_pair][0]]] + N_leaves_list[idx[pairs[min_pair][1]]]]))


    # List with all the previous min{d_ij}
    var_dij_history.append(dij_list[min_pair])

    linkage_list.append([idx[pairs[min_pair][0]], idx[pairs[min_pair][1]], min_tuple[0], N_leaves_list[-1]])


    logger.debug(f"------------------------------------------------------------")
    logger.debug(f"const_list= {const_list}")
    logger.debug(f"const_list[pairs[min_pair]]= {const_list[pairs[min_pair]]}")
    logger.debug(f"np.sum(const_list[pairs[min_pair]],axis=0) = {np.sum(const_list[pairs[min_pair]],axis=0)}")
    logger.debug(f"const_list[0] = {const_list[0]}")

    new_list = np.reshape(
      np.append(np.delete(const_list, pairs[min_pair], 0), [np.sum(const_list[pairs[min_pair]], axis=0)]), (-1, 2))
    logger.debug(f"New list =  {new_list}")

    jet_content = np.concatenate((jet_content, [np.sum(const_list[pairs[min_pair]], axis=0)]), axis=0)


    # Add a new key to the tree dictionary
    tree_dic[Nconst + Nparent] = idx[pairs[min_pair]]
    logger.debug(f"tree_dic = {tree_dic}")
    logger.debug(f"------------------------------------------------------------")

    # Delete the merged nodes
    idx = np.concatenate((np.delete(idx, pairs[min_pair]), [Nconst + Nparent]), axis=0)
    logger.debug(f"idx = {idx}")

    return new_list, var_dij_history, tree_dic, idx, jet_content, N_leaves_list, linkage_list






def _traverse(
        root,
        jet_nodes,
        tree_dic=None,
        Nleaves=None,
        dendrogram=True,
):
    """
    This function call the recursive function _traverse_rec to make the trees starting from the root
    :param root: root node id
    :param jet_nodes: array with the momentum of all the nodes of the jet tree (both leaves and inners).
    :param tree_dic: dictionary that has the node id of a parent as a key and a list with the id of the 2 children as the values
    :param Nleaves: Number of constituents (leaves)
    :param dendrogram: bool. If True, then return tree_ancestors list.

    :return:
    - tree: Reclustered tree structure.
    - content: Reclustered tree momentum vectors
    - node_id:   list where leaves idxs are added in the order they appear when we traverse the reclustered tree (each number indicates the node id
    that we picked when we did the reclustering.). However, the idx value specifies the order in which the leaf nodes appear when traversing the origianl jet (e.g. truth level) jet . The value here is an integer between 0 and Nleaves.
    So if we went from truth to kt algorithm, then in the truth tree the leaves go as [0,1,2,3,4,,...,Nleaves-1]
    - tree_ancestors: List with one entry for each leaf of the tree, where each entry lists all the ancestor node ids when traversing the tree from the root to the leaf node.

    """

    tree = []
    content = []
    node_id = []
    tree_ancestors = []

    _traverse_rec(
    root,
    -1,
    False,
    tree,
    content,
    jet_nodes,
    tree_dic=tree_dic,
    Nleaves=Nleaves,
    node_id=node_id,
    ancestors=[],
    tree_ancestors=tree_ancestors,
    dendrogram=dendrogram,
    )

    return tree, content, node_id, tree_ancestors






def _traverse_rec(
        root,
        parent_id,
        is_left,
        tree,
        content,
        jet_nodes,
        tree_dic=None,
        Nleaves=None,
        node_id=None,
        ancestors=None,
        tree_ancestors=[],
        dendrogram=False,
):
    """
    Recursive function to build the jet tree structure.
    :param root: parent node momentum
    :param parent_id: parent node idx
    :param is_left: bool.
    :param tree: List with the tree
    :param content: List with the momentum vectors
    :param jet_nodes: array with the momentum of all the nodes of the jet tree (both leaves and inners).
    :param tree_dic: dictionary that has the node id of a parent as a key and a list with the id of the 2 children as the values
    :param Nleaves: Number of constituents (leaves)
    :param node_id: list where leaves idxs are added in the order they appear when we traverse the reclustered tree (each number indicates the node id
    that we picked when we did the reclustering.). However, the idx value specifies the order in which the leaf nodes appear when traversing the truth level jet . The value here is an integer between 0 and Nleaves.
    So if we went from truth to kt algorithm, then in the truth tree the leaves go as [0,1,2,3,4,,...,Nleaves-1]
    :param ancestors: 1 entry of tree_ancestors (there is one for each leaf of the tree). It is appended to tree_ancestors.
    :param tree_ancestors: List with one entry for each leaf of the tree, where each entry lists all the ancestor node ids when traversing the tree from the root to the leaf node.
    :param dendrogram: bool. If True, append ancestors to tree_ancestors list.
    """


    id = len(tree) // 2
    if parent_id >= 0:
        if is_left:
            tree[2 * parent_id] = id  # We set the location of the lef child in the content array. So the left child momentum will be content[tree[2 * parent_id]]
        else:
            tree[2 * parent_id + 1] = id  # We set the location of the right child in the content array. So the right child will be content[tree[2 * parent_id+1]]
    """"(With each 4-vector we increase the content array by one element and the tree array by 2 elements. But then we take id=tree.size()//2, so the id increases by 1. The left and right children are added one after the other.)"""


    # Insert 2 new nodes to the vector that constitutes the tree.
    # In the next iteration we will replace this 2 values with the location of the parent of the new nodes
    tree.append(-1)
    tree.append(-1)

    # Fill the content vector with the values of the node
    content.append(jet_nodes[root])

    new_ancestors = None
    if dendrogram:
        new_ancestors = np.copy(ancestors)
        logger.debug(f" ancestors before = {ancestors}")
        new_ancestors = np.append(new_ancestors, root)  # Node ids in terms of the truth jet dictionary
        logger.debug(f" ancestors after = {ancestors}")


    # We move from the root down until we get to the leaves. We do this recursively
    if root >= Nleaves:

        children = tree_dic[root]
        logger.debug(f"Children = {children}")

        L_idx = children[0]
        R_idx = children[1]


        _traverse_rec(L_idx, id,
                      True,
                      tree,
                      content,
                      jet_nodes,
                      tree_dic,
                      Nleaves=Nleaves,
                      node_id=node_id,
                      ancestors=new_ancestors,
                      dendrogram=dendrogram,
                      tree_ancestors=tree_ancestors,
                      )

        _traverse_rec(R_idx,
                      id,
                      False,
                      tree,
                      content,
                      jet_nodes,
                      tree_dic,
                      Nleaves=Nleaves,
                      node_id=node_id,
                      ancestors=new_ancestors,
                      dendrogram=dendrogram,
                      tree_ancestors=tree_ancestors,
                      )

    # If not then its a leaf
    else:
        node_id.append(root)
        if dendrogram:
            tree_ancestors.append(new_ancestors)
            logger.debug(f"tree_ancestors= {tree_ancestors}")



if __name__== "__main__":

  input_dir = '../data/'
  input_jet = 'tree_0_truth'
  with open(input_dir + str(input_jet) + '.pkl', "rb") as fd:
    truth_jet = pickle.load(fd, encoding='latin-1')[0]

  jet_name = ('_').join(input_jet.split('_')[-3:-1])
  truth_jet["name"] = jet_name

  reclusterKt = recluster(truth_jet, alpha=1)