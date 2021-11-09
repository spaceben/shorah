import pytest
import filecmp
import os
import glob
from shorah import b2w, tiling
import math
import subprocess

p = os.path.dirname(__file__)

def _collect_files(base_path):
    spec_files = [os.path.basename(x) for x in glob.glob(os.path.join(base_path, '*.reads.fas'))]
    spec_files.extend(['coverage.txt', 'reads.fas'])
    return spec_files

@pytest.mark.parametrize("spec_dir,alignment_file,reference_file,region,window_length,overlap_factor,win_min_ext", [
    ("data_1", "test_aln.cram", "test_ref.fasta", "HXB2:2469-3713", 201, 3, 0.85), 
    ("data_1", "test_aln.cram", "test_ref.fasta", "HXB2:2469-3713", 204, 3, 0.85),
    ("data_1", "test_aln.cram", "test_ref.fasta", "HXB2:2469-3713", 200, 4, 0.85),  
    ("data_2", "REF_aln.bam", "cohort_consensus.fasta", "HXB2:2508-3676", 201, 3, 0.85),
    ("data_1", "test_aln.cram", "test_ref.fasta", "HXB2:2469-3713", 201, 3, 0.75), 
    ("data_1", "test_aln.cram", "test_ref.fasta", "HXB2:2469-3713", 200, 4, 0.65),
], indirect=["spec_dir"])
def test_cmp_raw(spec_dir, alignment_file, reference_file, region, window_length,overlap_factor, win_min_ext):
    assert window_length > 0 and window_length%overlap_factor == 0
    minimum_overlap = math.floor(window_length * win_min_ext)
    maximum_reads = math.floor(1e4 / window_length) # TODO why divide?
    incr = window_length//overlap_factor

    original = subprocess.run(
        f"b2w -w {window_length} -i {incr} -m {minimum_overlap} -x {maximum_reads} -c 0 {alignment_file} {reference_file} {region}", 
        shell=True, check=True, cwd=os.path.join(p, spec_dir)
    )
    assert original.returncode == 0

    strategy = tiling.EquispacedTilingStrategy(window_length, incr, True)

    b2w.build_windows(
        alignment_file = os.path.join(p, spec_dir, alignment_file),
        region = region,
        tiling_strategy = strategy,
        minimum_overlap = minimum_overlap, 
        maximum_reads = maximum_reads,
        minimum_reads = 0
    )

    spec_files = _collect_files(os.path.join(p, spec_dir))

    match, mismatch, errors = filecmp.cmpfiles(
        os.path.join(p, spec_dir),
        p,
        spec_files,
        shallow=False
    ) 
    print(match)
    print(mismatch)
    print(errors)

    assert len(mismatch) == len(errors) == 0


@pytest.fixture
def spec_dir(request):
    yield request.param

    # execute after test function, cleanup
    spec_files = _collect_files(os.path.join(p, request.param))
    created_files = _collect_files(p)

    for file in spec_files:
        os.remove(os.path.join(p, request.param, file))
    for file in created_files:
        os.remove(os.path.join(p, file))
    