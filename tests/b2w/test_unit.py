import filecmp
import os
import glob
from shorah import b2w

def test_cmp_raw():
    window_length = 201
    win_min_ext = 0.85 

    b2w.build_windows(
        alignment_file = "data/test_aln.cram",
        region = "HXB2:2469-3713",
        window_length = 201, 
        incr = window_length//3, 
        minimum_overlap = window_length * win_min_ext, 
        maximum_reads = 1e4 / window_length, # TODO why divide?
        minimum_reads = 0
    )

    
    p = os.path.dirname(__file__)

    spec_files = [os.path.basename(x) for x in glob.glob(os.path.join(p, 'data/*.fas'))]
    spec_files.extend(['coverage.txt'])

    match, mismatch, errors = filecmp.cmpfiles(
        os.path.join(p, 'data'),
        p,
        spec_files,
        shallow=False
    ) 
    print(match)
    print(mismatch)
    print(errors)

    assert len(mismatch) == len(errors) == 0
