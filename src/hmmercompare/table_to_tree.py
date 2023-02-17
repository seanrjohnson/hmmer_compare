"""Generates Trees from distance matrices
    
    Input is a table of bitscores from hmmer_compare (or any other program that produces scores in a tsv format with columns: query, reference, score).

    Ouput is a newick tree file.

    input is first converted to a dense matrix where rows and columns are the names of hmms and the values are the scores.

    Distances are then calculated with the formula:
        dist_matrix[i,j] = 1 - ( dist_matrix[i,j] / min(max(score_matrix[i,:],score_matrix[:,j])) )
"""

import argparse
from os import PathLike
import sys
from pathlib import Path
from hmmercompare import __version__, RawAndDefaultsFormatter
from typing import List, Union, Tuple
import warnings
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform
import numpy as np
from collections import OrderedDict

SCORE_COLUMN="TREE_DIST"

#see: https://stackoverflow.com/questions/28222179/save-dendrogram-to-newick-format
#see also: https://github.com/biopython/biopython/blob/master/Bio/Phylo/TreeConstruction.py
#see also: https://dendropy.org/library/phylogeneticdistance.html?highlight=upgma#dendropy.calculate.phylogeneticdistance.PhylogeneticDistanceMatrix.upgma_tree
# https://docs.scipy.org/doc/scipy/reference/cluster.hierarchy.html

def hmmer_compare_to_matrix(input_file: Union[str, PathLike]) -> Tuple[List[str], List[List[float]]]: 
    """Converts output from hmmer_compare to a distance matrix that can be used to generate a tree.
    
    Arguments:
        input_file {Union[str, PathLike]} -- the path to the input file
        output_file {Union[str, PathLike]} -- the path to the output file

    Returns:
        Tuple[List[str], List[List[float]]] -- a tuple containing the names and the distance matrix
    """
    input_file = Path(input_file)

    f = open(input_file, "r")

    f.readline() #skip header

    values = dict() # dict[(i,j)]:score
    indicies = OrderedDict() #dict[name]:index
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) != 3: # this allows us to parse files with alignments in them
            continue
        query, reference, score = parts[0], parts[1], float(parts[2])
        values[(query,reference)] = score
        if query not in indicies:
            indicies[query] = len(indicies)
        if reference not in indicies:
            indicies[reference] = len(indicies)
    matrix = np.zeros((len(indicies),len(indicies)))
    for (i,j), score in values.items():
        matrix[indicies[i],indicies[j]] = score
        matrix[indicies[j],indicies[i]] = score

    # calculate score dist
    matrix = 1 - (matrix/(np.minimum( matrix.max(axis=1)[:,np.newaxis], matrix.max(axis=0)[np.newaxis,:] ) ) )

    return list(indicies.keys()) , matrix

def clean_name_for_newick(name):
  """
    converts all: " " (space), ";" (semicolon), ":" (colon), "," (comma), "()" (parentheses), "'" (quote) characters to "_" in a string
  """

  bad_chars = " ;:,()'"
  chars = ["_" if x in bad_chars else x for x in name]
  out = "".join(chars)
  if out != name:
    warnings.warn(f"sequence name {name} contains characters that don't play well with newick, so for newick export, it has been changed to {out}")
  return out

def get_tree_edges(node, parent_dist, leaf_names, parent_name=None, edges=None, internal_nodes=None) -> dict: #adapted from: https://stackoverflow.com/questions/28222179/save-dendrogram-to-newick-format
    """
    Convert scipy.cluster.hierarchy.to_tree()-output to an edge list

    node: output of scipy.cluster.hierarchy.to_tree()
    parent_dist: output of scipy.cluster.hierarchy.to_tree().dist
    leaf_names: list of leaf names
    parent_name: the name of the parent node
    edges: leave empty, this variable will be populated by an empty dict, and is used in recursion.
    
    returns:
        dictionary of edges: [(start,end)]:dist
        internal_nodes: names of internal nodes
    """
    #TODO: make this not recursive

    if edges is None:
        edges = dict()
        internal_nodes = list()

    if node.is_leaf():
        name = leaf_names[node.id]
        edges[(parent_name,name)] = parent_dist - node.dist
    else:
        name = f"__INTERNAL_NODE_{node.id}__"
        internal_nodes.append(name)
        if parent_name is not None:
            edges[(parent_name,name)] = parent_dist - node.dist
        get_tree_edges(node.get_left(), node.dist, leaf_names, parent_name=name, edges=edges, internal_nodes=internal_nodes)
        get_tree_edges(node.get_right(), node.dist, leaf_names, parent_name=name, edges=edges, internal_nodes=internal_nodes)
        return edges, internal_nodes


def get_newick(node, parent_dist, leaf_names, newick='') -> str: #from: https://stackoverflow.com/questions/28222179/save-dendrogram-to-newick-format
    """ 
    Convert scipy.cluster.hierarchy.to_tree()-output to Newick format.

    node: output of scipy.cluster.hierarchy.to_tree()
    parent_dist: output of scipy.cluster.hierarchy.to_tree().dist
    leaf_names: list of leaf names
    newick: leave empty, this variable is used in recursion.
    
    returns:
        tree in Newick format
    """
    #TODO: make this not recursive

    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parent_dist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parent_dist - node.dist, newick)
        else:
            newick = ");"
        newick = get_newick(node.get_left(), node.dist, leaf_names, newick=newick)
        newick = get_newick(node.get_right(), node.dist, leaf_names, newick=",%s" % (newick))
        newick = "(%s" % (newick)
        return newick



def build_tree(scores_table_file, algorithm, newick_file):
    """ 
        Build a tree from a distance matrix.

        scores_table_file: path to the scores table file
        algorithm: the algorithm to use to build the tree. Currently only supports "upgma"
        newick_file: path to write the newick file to

    """
    row_names, matrix = hmmer_compare_to_matrix(scores_table_file)
    if algorithm == "upgma":
        
        linkage_matrix = hierarchy.average(squareform(matrix))
    
    tree = hierarchy.to_tree(linkage_matrix, False)

    if newick_file is not None:
        names = list(map(clean_name_for_newick, row_names))
        newick_string = get_newick(tree, tree.dist, names)
        with open(newick_file, "w") as outfile:
            print(newick_string, file=outfile)


def main(argv):
    parser = argparse.ArgumentParser(description=f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)
    
    parser.add_argument('-i', '--input', type=str, required=True,
                        help="Hmmer compare output file.")

    parser.add_argument('-o', "--output", type=str, required=True, default=None, #TODO: what other kinds of output might be useful?
                        help="Write a newick file of the tree.")

    parser.add_argument('--algorithm', type=str, required=False, default="upgma", choices={"upgma"},
                        help="Which clustering algorithm to use.")

    params = parser.parse_args(argv)


    input_file = params.input
    

    # Run
    # params.metric,
    build_tree(input_file, params.algorithm, params.output)


if __name__ == '__main__':
    main(sys.argv[1:])