"""
Sentetik (bilinçli olarak çift pik yerleştirilmiş) bir .ab1 dosyası
üretip, analiz algoritmasının bu pozisyonları doğru tespit ettiğini
doğrulayan basit testler.

Çalıştırmak için (pytest kurulu ise):
    pip install pytest
    pytest tests/

Ya da pytest olmadan doğrudan:
    python3 tests/test_analyzer.py
"""

import os
import sys
import shutil
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "examples"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from make_fake_ab1 import build_fake_ab1  # noqa: E402
from sanger_trace_analyzer import load_sanger_trace, analyze_trace  # noqa: E402


class TestSangerAnalyzer(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_mixed_positions_detected(self):
        sequence = "ACGGTTCAGCTAGGCATTAGCCA"
        mixed = {8: ("A", 0.40), 15: ("T", 0.60)}

        ab1_path = build_fake_ab1(
            os.path.join(self.tmpdir, "synthetic.ab1"), sequence, mixed
        )

        trace = load_sanger_trace(ab1_path)
        df = analyze_trace(trace)

        flagged_indices = set(df[df["flag"] == "mixed_candidate"]["index"])

        self.assertEqual(
            flagged_indices,
            set(mixed.keys()),
            f"Beklenen şüpheli pozisyonlar {set(mixed.keys())}, bulunan: {flagged_indices}",
        )

    def test_normal_positions_not_flagged(self):
        """Hiç mixed pozisyon eklenmemiş bir dizide hiçbir pozisyon
        'mixed_candidate' olarak işaretlenmemeli."""
        sequence = "ACGTACGTACGT"
        ab1_path = build_fake_ab1(
            os.path.join(self.tmpdir, "no_mixed.ab1"), sequence, {}
        )

        trace = load_sanger_trace(ab1_path)
        df = analyze_trace(trace)

        self.assertEqual((df["flag"] == "mixed_candidate").sum(), 0)

    def test_reads_correct_sequence(self):
        sequence = "ACGGTTCAGCTAGGCATTAGCCA"
        ab1_path = build_fake_ab1(
            os.path.join(self.tmpdir, "seq_check.ab1"), sequence, {}
        )

        trace = load_sanger_trace(ab1_path)
        self.assertEqual(trace.bases, sequence)


if __name__ == "__main__":
    unittest.main()
