"""Aligns and calculates alignment scores between hmmer3 profiles 

Kind of like hhsearch except much slower and for hmmer3 profiles instead of hhsuite profiles.

DISCLAIMER: this script is under development, I think the algorithm is correctly implemented, but I might change the output format and various other details.

Adapted from pseudocode in:
Steinegger, Martin, Markus Meier, Milot Mirdita, Harald Vöhringer, Stephan J. Haunsberger, and Johannes Söding. “HH-Suite3 for Fast Remote Homology Detection and Deep Protein Annotation.” BMC Bioinformatics 20, no. 1 (September 14, 2019): 473. https://doi.org/10.1186/s12859-019-3019-7.

and

Söding, Johannes. “Protein Homology Detection by HMM–HMM Comparison.” Bioinformatics 21, no. 7 (April 1, 2005): 951–60. https://doi.org/10.1093/bioinformatics/bti125.

"""
import argparse
from dataclasses import dataclass
import sys
import numpy as np
import pyhmmer
from typing import List, Iterable, TextIO, Tuple, Union, Dict
import os
from pathlib import Path
import heapq
from multiprocessing import Pool
from hmmercompare import __version__, RawAndDefaultsFormatter
import warnings




def read_hmms(hmm_files:Iterable[Union[str,os.PathLike]]) -> Dict[str, Dict[str,pyhmmer.plan7.HMM]]:
    """
        hmm_files: a list of paths to .hmm files

        returns:
            a dict of dicts of pyhmmer HMM objects
                db_name: hmm_name: HMM

    """
    out = dict()
    for file in hmm_files:
        name = os.path.basename(Path(file).stem)
        
        hmmer_models = dict() 
        for model in pyhmmer.plan7.HMMFile(file):
            model_name = model.name.decode()
            if model_name in hmmer_models:
                warnings.warn(f"multiple hmms with the same name ({model_name}) in file: {file}, only one will be used.")
            hmmer_models[model_name] = model
        if name in out:
            raise RuntimeError(f"Multiple hmm files with the same name, please combine the hmms into a single file, or rename one of the files. This is important to avoid searching the same domains twice, and for naming the source databases.")
        out[name] = hmmer_models
    return out


#from Algorithm2 pseudocode in:
# Steinegger, Martin, Markus Meier, Milot Mirdita, Harald Vöhringer, Stephan J. Haunsberger, and Johannes Söding. “HH-Suite3 for Fast Remote Homology Detection and Deep Protein Annotation.” BMC Bioinformatics 20, no. 1 (September 14, 2019): 473. https://doi.org/10.1186/s12859-019-3019-7.

#TODO: cythonize

_SMM = 0
_SMI = 1
_SIM = 2
_SDG = 3
_SGD = 4
_TRACEBACK_DTYPE=np.uint64
_STOP_FLAG = np.iinfo(_TRACEBACK_DTYPE).max


def max2(sMM, sXY, layer1, layer2) -> Tuple[float, int]:
    if sMM > sXY:
        score = sMM
        bt = layer1
    else:
        score = sXY
        bt = layer2
    return score, bt


def max6(sMM: float, sMI: float, sIM: float, sDG: float, sGD: float, sSTOP: float = 0.0) -> Tuple[float, int]:
    """

    Args:
        sMM (float): 
        sMI (float): 
        sIM (float): 
        sDG (float): 
        sGD (float):         

        sSTOP (float, optional): Minimum possible alignment score. 0 for local alignment, -inf for global alignment. Defaults to 0. 

    Returns:
        Tuple[float, int]: score, backtrace_layer

    """
    score = float("-inf")
    bt = None
    
    if sSTOP > sMM:
        score = sSTOP
        bt = _STOP_FLAG
    else:
        score = sMM
        bt = _SMM
    
    if sIM > score:
        score = sIM
        bt = _SIM
    if sMI > score:
        score = sMI
        bt = _SMI
    if sDG > score:
        score = sDG
        bt = _SDG
    if sGD > score:
        score = sGD
        bt = _SGD
    
    return  score, bt
    

def Saa(q:np.array, r:np.array, background:np.array) -> float: #TODO: test
    """calculate similarity scores for two amino acid distributions

    Args:
        q (np.array): array size 20, amino acid probabilities
        r (np.array): array size 20, amino acid probabilities
        background (np.array): array size 20,amino acid background probabilities

    Returns:
        float: score
    """

    return np.log(np.sum(np.divide(np.multiply(q, r), background)))


def compare_hmmer(qhmm:pyhmmer.plan7.HMM, rhmm:pyhmmer.plan7.HMM) -> Tuple[float, np.array, Tuple[int,int,int]]:
    """Run the viterbi alignment algorithm to compare two plan7 HMM profiles

    Args:
        qhmm (pyhmmer.plan7.HMM): query HMM object
        rhmm (pyhmmer.plan7.HMM): reference (target) HMM object

    Raises:
        ValueError: if hmms are not of the same alphabet

    Returns:
        Tuple[float, np.array, Tuple[int,int,int]]: score, backtrace, max_index (query_idx, ref_idx, layer)
    """
    if qhmm.alphabet != rhmm.alphabet:
        raise ValueError(f"Error, cannot compare hmms with different alphabets: {qhmm.alphabet}, {rhmm.alphabet}.")
    # alphabet_K = qhmm.alphabet.K
    background = np.asarray(pyhmmer.plan7.Background(qhmm.alphabet).residue_frequencies) #TODO: could just store this as a constant

    # raw probabilities
    query_match_emissions=np.asarray(qhmm.match_emissions) # M + 1, K. Row 0 is entry probabilities, so emissions are unused, but it still needs to sum to 1.0, so it is [1.0, 0.0, 0.0, ...], but don't use it for any calcs
    reference_match_emissions=np.asarray(rhmm.match_emissions) # M + 1, K. Row 0 is entry probabilities, so emissions are unused, but it still needs to sum to 1.0, so it is [1.0, 0.0, 0.0, ...], but don't use it for any calcs
    
    # query_insert_emissions=np.asarray(qhmm.insert_emissions) # row 0 is insert emissions for left-gaps (before the first match emission) 
    # reference_insert_emissions=np.asarray(rhmm.insert_emissions) # row 0 is insert emissions for left-gaps (before the first match emission)

    with np.errstate(divide='ignore'):
        query_transitions=-1 * np.log(np.asarray(qhmm.transition_probabilities))  # M + 1, 7 probabilities are negative natural logs
        reference_transitions=-1 * np.log(np.asarray(rhmm.transition_probabilities)) # M + 1, 7 probabilities are negative natural logs
    
    # print(np.asarray(qhmm.transition_probabilities))
    TMM=0
    TMI=1
    TMD=2
    TIM=3
    TII=4
    TDM=5
    TDD=6
    # 0: Mn -> Mn+1 
    # 1: Mn -> In+1 
    # 2: Mn -> Dn+1 
    # 3: In -> Mn+1 
    # 4: In -> In+1 
    # 5: Dn -> Mn+1 
    # 6: Dn -> Dn+1 

    scores = np.zeros((query_match_emissions.shape[0], reference_match_emissions.shape[0],5),dtype=np.double) # dim3 = SMM, SMI, SIM, SDG, SGD, for SXY, X is query, Y is reference
    backtrace = np.zeros((query_match_emissions.shape[0], reference_match_emissions.shape[0],5),dtype=_TRACEBACK_DTYPE) #score_row, score_column, state_type (dim3 of scores array) TODO: should we account for alternative alignments, like by ties in the max fuctions?. TODO: what order of dimensions is best for caching?
    match_scores = np.zeros((query_match_emissions.shape[0], reference_match_emissions.shape[0]),dtype=np.double) # scores for match states only, for backtrace
    # values are flags to the previous state type. There is only one direction from each previous state, so directionality is implied in the state type

    backtrace[:,0,:] = _STOP_FLAG # stop backtrace at the left column
    backtrace[0,:,:] = _STOP_FLAG # stop backtrace at the top row
    # SMM = 0
    # SMI = 1
    # SIM = 2
    # SDG = 3
    # SGD = 4
    
    for qi in range(1,query_match_emissions.shape[0]):
        for ri in range(1,reference_match_emissions.shape[0]):

            # SMM
            scores[qi,ri,_SMM], tb_layer = \
                max6( # find the best previous state
                    scores[qi-1,ri-1,_SMM] - query_transitions[qi-1, TMM] - reference_transitions[ri-1, TMM], # subtracting negative logs is the same as adding logs, which is the same as multiplying unlogs
                    scores[qi-1,ri-1,_SMI] - query_transitions[qi-1, TMM] - reference_transitions[ri-1, TIM],
                    scores[qi-1,ri-1,_SIM] - query_transitions[qi-1, TIM] - reference_transitions[ri-1, TMM], 
                    scores[qi-1,ri-1,_SDG] - query_transitions[qi-1, TDM] - reference_transitions[ri-1, TMM], 
                    scores[qi-1,ri-1,_SGD] - query_transitions[qi-1, TMM] - reference_transitions[ri-1, TDM],
                 ) 
            backtrace[qi,ri,_SMM] = tb_layer
            #Saa_score = Saa(query_match_emissions[qi,:], reference_match_emissions[ri,:], background)
            #print(Saa_score)
            match_score = Saa(query_match_emissions[qi,:], reference_match_emissions[ri,:], background)
            match_scores[qi,ri] = match_score
            scores[qi,ri,_SMM] += match_score
            # could add a secondary structure term here if available
            
            # SGD
            scores[qi,ri,_SGD], tb_layer = \
                max2(
                    scores[qi,ri-1,_SMM] - reference_transitions[ri-1, TMD],
                    scores[qi,ri-1,_SGD] - reference_transitions[ri-1, TDD],
                    _SMM,
                    _SGD
                )
            backtrace[qi,ri,_SGD] = tb_layer

            # SIM
            scores[qi,ri,_SIM], tb_layer = \
                max2(
                    scores[qi,ri-1,_SMM] - query_transitions[qi, TMI] - reference_transitions[ri-1, TMM],
                    scores[qi,ri-1,_SIM] - query_transitions[qi, TII] - reference_transitions[ri-1, TMM],
                    _SMM,
                    _SIM
                )
            backtrace[qi,ri,_SIM] = tb_layer


            # SDG
            scores[qi,ri,_SDG], tb_layer = \
                max2(
                    scores[qi-1,ri,_SMM] - query_transitions[qi-1, TMD],
                    scores[qi-1,ri,_SDG] - query_transitions[qi-1, TDD],
                    _SMM,
                    _SDG
                )
            backtrace[qi,ri,_SDG] = tb_layer

            # SMI
            scores[qi,ri,_SMI], tb_layer = \
                max2(
                    scores[qi-1,ri,_SMM] - query_transitions[qi-1, TMM] - reference_transitions[ri, TMI],
                    scores[qi-1,ri,_SMI] - query_transitions[qi-1, TMM] - reference_transitions[ri, TII],
                    _SMM,
                    _SMI
                )
            backtrace[qi,ri,_SMI] = tb_layer


    max_index = np.unravel_index(np.argmax(scores),scores.shape)  
    score = scores[max_index]
    

    return score, backtrace, max_index, match_scores

def traceback(qhmm,rhmm,backtrace,trace_start,match_scores):
    POSITION_PADDING = 10
    
    #alphabet_symbols = qhmm.alphabet.symbols # ACDEFGHIKLMNPQRSTVWY-BJZOUX*~
    qcons = qhmm.consensus
    rcons = rhmm.consensus
    qend = trace_start[0] # int 
    rend = trace_start[1] # int
    qstart = qend
    rstart = rend

    ptr = trace_start # 3-tuple q_pos, r_pos, tb_layer

    prev = backtrace[ptr]
    qline=list()
    midline=list()
    rline=list()
    # SMM = 0
    # SMI = 1
    # SIM = 2
    # SDG = 3
    # SGD = 4

    # hhsearch style alignment output symbols TODO: are these cutoffs calibrated the same way with hmmer3, or do I need to adjust them somehow?
    # https://github.com/soedinglab/hh-suite/wiki#hmm-hmm-pairwise-alignments
    # = : column score below -1.5
    # - : column score between -1.5 and -0.5
    # . : column score between -0.5 and +0.5
    # + : column score between +0.5 and +1.5
    # | : column score above   +1.5


    while True:
        ptr_layer = ptr[2]
        qpos = ptr[0]-1 #-1 because the backtrace array has an extra row at the front
        rpos = ptr[1]-1 #-1 because the backtrace array has an extra column at the front
        if ptr_layer == _SMM:
            qchar = qcons[qpos]
            rchar = rcons[rpos]
            qline.append(qchar)
            rline.append(rchar)
            next_ptr = (ptr[0]-1, ptr[1]-1, prev)
            qstart = ptr[0]
            rstart = ptr[1]

            score = match_scores[qpos+1,rpos+1]
            if score < -1.5:
                midline.append("=")
            elif score < -0.5:
                midline.append("-")
            elif score < 0.5:
                midline.append(".")
            elif score < 1.5:
                midline.append("+")
            else:
                midline.append("|")
        elif ptr_layer == _SMI:
            qline.append(qcons[qpos])
            rline.append("-")
            midline.append(" ")
            qstart = ptr[0]
            next_ptr = (ptr[0]-1, ptr[1], prev)
        elif ptr_layer == _SIM:
            qline.append("-")
            rline.append(rcons[rpos])
            midline.append(" ")
            rstart = ptr[1]
            next_ptr = (ptr[0], ptr[1]-1, prev)
        elif ptr_layer == _SDG:
            qline.append(qcons[qpos])
            rline.append("-")
            midline.append(" ")
            qstart = ptr[0]
            next_ptr = (ptr[0]-1, ptr[1], prev)
        elif ptr_layer == _SGD:
            qline.append("-")
            rline.append(rcons[rpos])
            midline.append(" ")
            rstart = ptr[1]
            next_ptr = (ptr[0], ptr[1]-1, prev)
        
        if prev == _STOP_FLAG:
            break
        ptr = next_ptr
        prev = backtrace[ptr]
    return str(qstart).ljust(POSITION_PADDING) + ''.join(reversed(qline)) + str(qend).rjust(POSITION_PADDING) + "\n" + \
           ' '*POSITION_PADDING + ''.join(reversed(midline)) + ' ' * POSITION_PADDING + "\n" + \
           str(rstart).ljust(POSITION_PADDING) + ''.join(reversed(rline)) + str(rend).rjust(POSITION_PADDING)

@dataclass
class HmmerHit():
    score: float
    query_name: str
    reference_name: str
    alignment: str
    
    def __eq__(self, other):
        return self.score == other.score

    def __ne__(self, other):
        return self.score != other.score

    def __lt__(self, other):
        return self.score < other.score

    def __le__(self, other):
        return self.score <= other.score

    def __gt__(self, other):
        return self.score > other.score

    def __ge__(self, other):
        return self.score >= other.score


class _hmmer_compare_worker():
    def __init__(self, hmmer_targets, alignment=False, k=None, score_cutoff=float("-inf")):
        self.k = k
        self.hmmer_targets = hmmer_targets
        self.alignment = alignment
        self.score_cutoff = score_cutoff
    
    def __call__(self, input_profile):
        out_heap = []
        for target_dataset in self.hmmer_targets.values():
            for target_profile in target_dataset.values():
                score, backtrace, max_index, match_scores = compare_hmmer(input_profile, target_profile)
                if score >= self.score_cutoff:
                    if self.alignment:
                        alignment = traceback(input_profile,target_profile,backtrace,max_index, match_scores)
                    else:
                        alignment = None
                    
                    result = HmmerHit(round(score,2), input_profile.name.decode(), target_profile.name.decode(), alignment)
                    if (self.k is None) or (len(out_heap) < self.k):
                        heapq.heappush(out_heap, result)
                    else:
                        heapq.heappushpop(out_heap, result)
        out_heap.sort(reverse=True)
        return out_heap



def compare_hmmer_files(query_files:Iterable[str], reference_files:Iterable[str], out_handle:TextIO, score_cutoff:float, alignments:bool, k:int, cpu:int):
    references = read_hmms(hmm_files=reference_files) # list of lists of pyhmmer hmm objects

    worker = _hmmer_compare_worker(references, alignments, k, score_cutoff)

    sep="\t"
    print(sep.join(("query","reference","score")), file=out_handle) #TODO: how to write the alignment?
    
    for file in query_files:
        # file_name = os.path.basename(Path(file).stem)
        with Pool(processes=cpu) as pool:
            for hits in pool.imap_unordered(worker, pyhmmer.plan7.HMMFile(file)):
                for hit in hits:
                    print(sep.join( (hit.query_name,hit.reference_name,f"{round(hit.score,2)}") ), file=out_handle)
                    if alignments:
                        print(hit.alignment, file=out_handle)
                        print("\n\n", file=out_handle)


def main(argv):
    parser = argparse.ArgumentParser(f"\nversion: {__version__}\n\n" + __doc__, formatter_class=RawAndDefaultsFormatter)
    parser.add_argument('-i', '--input', type=str, required=True, nargs='+',
                        help="Input query files. One or more hmm text files with one or more hmmer3 profiles.")
    parser.add_argument('-r', "--reference", type=str, required=True, nargs='+',
                        help="Reference files. One or more hmm text files with one or more hmmer3 profiles.")

    parser.add_argument('--score_cutoff', type=float, default = 0,
                        help="Report alignments with scores greater than or equal to this.") #TODO: what is a reasonable cutoff?

    parser.add_argument('-k', type=int, required=False, default=None,
                        help="Include at most this many of the top hits for each query. Default: Include all hits.")

    parser.add_argument('-o', '--output', type=str, default=None,
                        help="File to write the scores and alignments to.")

    parser.add_argument('--alignments', action='store_true', default=False,
                        help="when activated, will write the alignments to the output.")
    
    parser.add_argument('--cpu', type=int, default=8, required=False,
                        help="how many cpu threads to use. Default: 8")


    params = parser.parse_args(argv)

    if params.output is None:
        out = sys.stdout
    else:
        out = open(params.output, "w")

    compare_hmmer_files(params.input, params.reference, out, params.score_cutoff, params.alignments, params.k, params.cpu)

    if params.output is not None:
        out.close()

if __name__ == '__main__':
    main(sys.argv[1:])
