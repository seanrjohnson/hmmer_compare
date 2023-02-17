import tempfile
from hmmercompare import hmmer_compare

#TODO: better tests!

def test_hmmer_compare_1(shared_datadir):

    with tempfile.TemporaryDirectory() as output_dir:
        #output_dir = "tmp_out"
        out_path = output_dir + f"/out_scores.tsv"
        hmmer_compare.main(["-i", str(shared_datadir / "pdonr_hmms.hmm"),  "-r", str(shared_datadir / "pdonr_hmms.hmm"), "-o", out_path, "--alignment", "--score_cutoff", "13"])

