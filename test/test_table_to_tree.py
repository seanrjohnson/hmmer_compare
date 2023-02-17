import tempfile
from hmmercompare import table_to_tree
from helpers import compare_files

#TODO: better tests!

def test_hmmer_compare_1(shared_datadir):

    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "tmp_out"
        out_path = output_dir + f"/out_scores.tsv"
        table_to_tree.main(["-i", str(shared_datadir / "tree_scores.tsv"), "-o", out_path])
        compare_files(out_path, shared_datadir / "pdonr.newick")
