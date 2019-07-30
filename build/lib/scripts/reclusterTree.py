import sys
# sys.path.append("..")
# from IPython.display import display

import numpy as np
import logging
import pickle
import itertools

# from sklearn.metrics import roc_auc_score
# from sklearn.preprocessing import RobustScaler
# from sklearn.utils import check_random_state

# %matplotlib inline
# import matplotlib.pyplot as plt
# plt.rcParams["figure.figsize"] = (8,8)


def recluster(input_jet, alpha=None):
  """
  Uses helper functions to get the leaves of a jet, recluster them following some algorithm, create the new tree for the
  algorithm, make a jet dictionary and save it.

   Create a dictionary with all the jet tree info (topology, constituents features: eta, phi, pT, E)
  - tree: array witht the clustering history
  - content: array with the pseudojets momentum

  Args:
  - input_jet: any jet dictionary with the clustering history.
  - alpha: defines the clustering algorithm. alpha={-1,0,1} defines the {anti-kt, CA and kt} algorithms respectively.
  """

  def _rec(jet, parent, node_id, outers_list):
    """
    Recursive function to get a list of the leaves
    """
    if jet["tree"][node_id, 0] == -1:
      outers_list.append(jet["content"][node_id])
    else:
      _rec(jet, node_id, jet["tree"][node_id, 0], outers_list)
      _rec(jet, node_id, jet["tree"][node_id, 1], outers_list)

    return outers_list

  outers = []
  jet_const = np.asarray(_rec(input_jet, -1, input_jet["root_id"], outers))

  raw_tree, \
  idx, \
  jet_content, \
  root_node, \
  Nconst, \
  N_leaves_list, \
  linkage_list = ktAntiktCA(jet_const, alpha=alpha)

  tree, \
  content, \
  node_id, \
  tree_ancestors = _traverse(root_node,
                             jet_content,
                             tree_dic=raw_tree,
                             root_idx=None,
                             Nleaves=Nconst,
                             )

  jet = {}
  jet["root_id"] = 0
  jet["tree"] = np.asarray(tree).reshape(-1, 2)
  jet["content"] = np.asarray([np.asarray(c) for c in content]).reshape(-1, 2)

  #     print(jet)

  # Save reclustered tree
  out_dir = 'data/'
  algo = str(jet_dic["name"]) + '_' + str(alpha)
  out_filename = out_dir + str(algo) + '.pkl'
  print('out_filename=', out_filename)
  with open(out_filename, "wb") as f:
    pickle.dump(jet, f, protocol=2)

  return node_id, linkage_list, Nconst, tree_ancestors

def ktAntiktCA(const_list, alpha=None):
  """
  Runs the dijMinPair function level by level until we reach the root of the tree.

  Args:
      - const_list: jet constituents (i.e. the leaves of the tree)
      - alpha: defines the clustering algorithm. alpha={-1,0,1} defines the {anti-kt, CA and kt} algorithms respectively.

  Returns:
      - const_list: nodes list after deleting the constituents that are merged and adding the new pseudojet in all levels.
        So this should only have the root of the tree.
      - tree_dic: dictionary that has the node id of a parent as a key and a list with the id of the 2 children as the values
      - idx: array that stores the node id
       (the node id determines the location of the momentum vector of a pseudojet in the jet_content array)
        of the pseudojets that are in the current const_list array. It has the same elements as the const_list (they get updated
        level by level).
      - jet_content: array with the momentum of all the nodes of the jet tree (both leaves and inners).
      - root_node: root node id
      - Nconst: Number of leaves of the jet
  """
  Nconst = len(const_list)

  root_node = 2 * Nconst - 2
  #     print('Root node= (N constituents + N parent) =', root_node)

  idx = np.arange(Nconst)

  # List that given a node idx, stores for that idx, the number of leaves for the branch below that node.
  # It is initialized only with the tree leaves
  N_leaves_list = np.ones((Nconst))

  # List that will have the info to do the dendrogram heat map plots
  linkage_list = []

  dij_hist = []
  tree_dic = {}

  const_list = np.asarray(const_list)
  #     const_list = np.reshape(const_list,len(const_list))
  jet_content = const_list

  #     print('const_list =', const_list)
  #     print('CONST LIST SUM = ',np.sum(const_list))
  #     print('---'*20)

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
    -Calculate all d_ij and get the minimum.
    -Update the constituents list by deleting the constituents that are merged and adding the new pseudojet

    Note: We refer to both leaves and inner nodes as pseudojets.

    Args:
        - const_list: constituents list for the current level (i.e. deleting the constituents that are merged and
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
        - linkage_list: A (Nconst-1) by 4 matrix that is used as input to build the dendrogram heat map plots.
           At the -th iteration, clusters with indices Z[i, 0] and Z[i, 1] are combined to form cluster .
           A cluster with an index less than  corresponds to one of the  original observations.
           The distance between clusters Z[i, 0] and Z[i, 1] is given by Z[i, 2].
           The fourth value Z[i, 3] represents the number of original observations in the newly formed cluster.

    Returns:
        -new_list: new const_list after deleting the constituents that are merged and adding the new pseudojet in the current level.
        -var_dij_history
        - tree_dic
        - idx
        - jet_content

    """

    # Get all possible pairings
    pairs = np.asarray(list(itertools.combinations(np.arange(len(const_list)), 2)))

    #     print('Const list all=', const_list)
    const_list_pt = np.absolute([element[0] for element in const_list])
    #     print('const_list_pt =', const_list_pt)
    #     print('const_list_pt[pairs] =',const_list_pt[pairs])

    # Get all dij at each level: dij=min(pTi^(2\alpha),pTj^(2\alpha)) * [arccos((pi.pj)/|pi|*|pj|)]^2
    dij_list = [(np.sort((const_list_pt[pairs][k]) ** (2 * alpha))[0] * \
                 (np.arccos(np.sum(const_list[pairs][k][0] * const_list[pairs][k][1]) /
                            (np.sqrt(np.sum(const_list[pairs][k][0] ** 2)) * np.sqrt(
                              np.sum(const_list[pairs][k][1] ** 2))))) ** 2, k) for k in range(len(const_list[pairs]))]
    #     print('dij_list =',dij_list)

    # Get pair index (in pairs list) with min dij
    min_tuple = sorted(dij_list, key=lambda x: x[0])[0]
    min_pair = min_tuple[1]
    #     print('min_pair=',pairs[min_pair])

    # List that given a node idx, stores for that idx, the number of leaves for the branch below that node.
    N_leaves_list = np.concatenate(
      (N_leaves_list, [N_leaves_list[idx[pairs[min_pair][0]]] + N_leaves_list[idx[pairs[min_pair][1]]]]))
    #     print('idx[pairs[min_pair]=', idx[pairs[min_pair]])

    # List with all the previous min{d_ij}
    var_dij_history.append(dij_list[min_pair])

    linkage_list.append([idx[pairs[min_pair][0]], idx[pairs[min_pair][1]], min_tuple[0], N_leaves_list[-1]])

    #     print('+++++'*3)
    #     print('const_list=', const_list)
    #     print('const_list[pairs[min_pair]]=', const_list[pairs[min_pair]])
    #     print('np.sum(const_list[pairs[min_pair]],axis=0) =', np.sum(const_list[pairs[min_pair]],axis=0))
    #     print('const_list[0]', const_list[0])

    new_list = np.reshape(
      np.append(np.delete(const_list, pairs[min_pair], 0), [np.sum(const_list[pairs[min_pair]], axis=0)]), (-1, 2))
    #     print('New list =', new_list)

    jet_content = np.concatenate((jet_content, [np.sum(const_list[pairs[min_pair]], axis=0)]), axis=0)

    # Add a new key to the tree dictionary
    tree_dic[Nconst + Nparent] = idx[pairs[min_pair]]
    #     print('tree_dic = ', tree_dic)
    #     print('---'*5)

    # Delete the merged nodes
    idx = np.concatenate((np.delete(idx, pairs[min_pair]), [Nconst + Nparent]), axis=0)
    #     idx.append(Nconst+Nparent)
    #     print('idx =',idx)

    return new_list, var_dij_history, tree_dic, idx, jet_content, N_leaves_list, linkage_list


# -------------------------------------------------------------------------------------------------------------
# This function call the recursive function to make the trees starting from the root
def _traverse(
        root,
        jet_nodes,
        tree_dic=None,
        root_idx=None,
        Nleaves=None,
        dendrogram=True,
):  # root should be a fj.PseudoJet

  tree = []
  content = []
  node_id = []
  #     ancestors=[]
  tree_ancestors = []

  #   sum_abs_charge=0
  _traverse_rec(
    root,
    -1,
    False,
    tree,
    content,
    jet_nodes,
    tree_dic=tree_dic,
    root_idx=root_idx,
    Nleaves=Nleaves,
    node_id=node_id,
    ancestors=[],
    tree_ancestors=tree_ancestors,
    dendrogram=dendrogram,
  )  # We start from the root=jet 4-vector

  #     print('-- tree_ancestors=',tree_ancestors)
  return tree, content, node_id, tree_ancestors

# -------------------------------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------------------------------
# Recursive function to access fastjet clustering history and make the tree. We will call this function below in _traverse.
def _traverse_rec(
        root,
        parent_id,
        is_left,
        tree,
        content,
        jet_nodes,
        tree_dic=None,
        root_idx=None,
        Nleaves=None,
        node_id=None,
        ancestors=None,
        tree_ancestors=[],
        dendrogram=False,
):  # root should be a fj.PseudoJet

  '''
  Args:
  root: index of a node
  '''

  id = len(tree) // 2
  #     print('id=', id)
  if parent_id >= 0:
    if is_left:
      tree[
        2 * parent_id] = id  # We set the location of the lef child in the content array of the 4-vector stored in content[parent_id]. So the left child will be content[tree[2 * parent_id]]
    else:
      tree[
        2 * parent_id + 1] = id  # We set the location of the right child in the content array of the 4-vector stored in content[parent_id]. So the right child will be content[tree[2 * parent_id+1]]
  #  This is correct because with each 4-vector we increase the content array by one element and the tree array by 2 elements. But then we take id=tree.size()//2, so the id increases by 1. The left and right children are added one after the other.

  # -------------------------------
  # We insert 2 new nodes to the vector that constitutes the tree. In the next iteration we will replace this 2 values with the location of the parent of the new nodes
  tree.append(-1)
  tree.append(-1)

  #     We fill the content vector with the values of the node
  content.append(jet_nodes[root])

  new_ancestors = None
  if dendrogram:
    new_ancestors = np.copy(ancestors)
    #         print(' ancestors before=', ancestors)
    #         new_ancestors = np.append(new_ancestors,id) # Node ids in terms of the reclustered jet dictionary
    new_ancestors = np.append(new_ancestors, root)  # Node ids in terms of the truth jet dictionary
  #         print('ancestors=', ancestors)
  #         print('New ancestors =', new_ancestors)
  #         print('root=', id)
  #         print('Nleaves= ', Nleaves)
  #         print('---'*5)

  #   content.append(root.py())
  #   content.append(root.pz())
  #   content.append(root.e())

  # --------------------------------------
  # We move from the root down until we get to the leaves. We do this recursively

  #     If not then its a leaf
  if root >= Nleaves:

    children = tree_dic[root]
    #         print('Children = ', children)

    L_idx = children[0]
    R_idx = children[1]

    # ------------------------------
    # Call the function recursively

    _traverse_rec(L_idx, id,
                  True,
                  tree,
                  content,
                  jet_nodes,
                  tree_dic,
                  root_idx=L_idx,
                  Nleaves=Nleaves,
                  node_id=node_id,
                  ancestors=new_ancestors,
                  dendrogram=dendrogram,
                  tree_ancestors=tree_ancestors,
                  )  # pieces[0] is the left child

    _traverse_rec(R_idx,
                  id,
                  False,
                  tree,
                  content,
                  jet_nodes,
                  tree_dic,
                  root_idx=R_idx,
                  Nleaves=Nleaves,
                  node_id=node_id,
                  ancestors=new_ancestors,
                  dendrogram=dendrogram,
                  tree_ancestors=tree_ancestors,
                  )  # pieces[1] is the right child


  else:
    node_id.append(
      root)  # Fill idx list with the index that specifies the order in which the leaf nodes appear when traversing
    #  the truth level jet . The root value here is an integer between 0 and Nleaves.
    # And these numbers are added in the order that they appear when we traverse the tree. Each number indicates the node id
    # that we pick when we did the reclustering.
    # So if we went from truth to kt algorithm, then in the truth tree the leaves go as [0,1,2,3,4,,...,Nleaves-1]
    if dendrogram:
      tree_ancestors.append(new_ancestors)
#             print('tree_ancestors=', tree_ancestors)