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
    seen_mutations = set()
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
                seen_mutations.add(key)
            else:
                assert row[key] == "0", f"Erroneously found {key} in {sample}: {row}"
        assert f"{sample} found" in comment, f"Missing '{sample} found' from {comment!r}"
        assert sample_found, f"Failed to find {sample} in {row}"

    expected_mutations = set(reader._fieldnames)
    expected_mutations.remove("sample")
    expected_mutations.remove("comment")

    missed_mutatations = expected_mutations.difference(seen_mutations)
    assert len(missed_mutatations) == 0, f"No tests for mutations: {missed_mutatations}"


if __name__ == "__main__":
    tmpdir = pathlib.Path(tempfile.mkdtemp())
    try:
        integration_run(tmpdir)
    finally:
        shutil.rmtree(tmpdir)
