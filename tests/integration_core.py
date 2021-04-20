"""Integration tests"""

import csv
import io
import pathlib
import shutil
import subprocess
import tempfile

def integration_run(tmp_path):
    datadir = pathlib.Path(__file__).parent / "data"
    cmd = [
        "covid-spike-classification",
        "--input-format", "fasta",
        "--outdir", tmp_path,
        "--stdout",
        "--silence-warnings",
        "--quiet",
        datadir,
    ]

    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    outs, _ = proc.communicate(timeout=30)

    assert len(outs) > 0, "command generated no output"
    csv_handle = io.StringIO(outs.decode("utf-8"))
    reader = csv.DictReader(csv_handle)
    for row in reader:
        sample = row['sample']
        comment = row['comment']
        sample_found = False
        for key in row.keys():
            if key in ("sample", "comment"):
                continue
            if key == sample:
                assert row[key] == "1", f"Failed to find {sample}"
                sample_found = True
            else:
                assert row[key] == "0", f"Erroneously found {key} in {sample}: {row}"
        assert f"{sample} found" in comment, f"Missing '{sample} found' from {comment!r}"
        assert sample_found, f"Failed to find {sample} in {row}"


if __name__ == "__main__":
    tmpdir = pathlib.Path(tempfile.mkdtemp())
    try:
        integration_run(tmpdir)
    finally:
        shutil.rmtree(tmpdir)
