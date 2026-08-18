import gzip
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "fastq_run_qc.py"
SPEC = importlib.util.spec_from_file_location("fastq_run_qc", SCRIPT)
QC = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = QC
SPEC.loader.exec_module(QC)


class FastqRunQcTests(unittest.TestCase):
    def write_fastq(self, path: Path, records: list[tuple[str, str, str]]) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        handle_context = gzip.open(path, "wt") if path.name.endswith(".gz") else path.open("wt")
        with handle_context as handle:
            for header, sequence, quality in records:
                handle.write(f"@{header}\n{sequence}\n+\n{quality}\n")

    def test_directory_qc_and_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self.write_fastq(
                root / "barcode01" / "chunk.fastq",
                [("read1", "AAAA", "IIII"), ("read2", "ACGTACGT", "+" * 8)],
            )
            self.write_fastq(
                root / "barcode02" / "chunk.fastq.gz",
                [("read3", "A" * 10, "5" * 10)],
            )

            inputs = QC.discover_fastqs(root)
            self.assertEqual([sample for _, sample in inputs], ["barcode01", "barcode02"])
            stats = QC.analyse_fastqs(inputs, "path", 0, 0)
            barcode01 = QC.metrics_for(stats["barcode01"], 2, 0.5)
            self.assertEqual(barcode01["reads"], 2)
            self.assertEqual(barcode01["bases"], 12)
            self.assertEqual(barcode01["median_length"], 6)
            self.assertEqual(barcode01["n50"], 8)
            self.assertAlmostEqual(barcode01["q10_percent"], 100.0)
            self.assertAlmostEqual(barcode01["q15_percent"], 50.0)

            output = root / "results"
            self.assertEqual(
                QC.main([str(root), "--barcode-source", "path", "--output-dir", str(output), "--progress-every", "0"]),
                0,
            )
            self.assertTrue((output / "run_qc_dashboard.html").is_file())
            self.assertTrue((output / "run_qc_summary.tsv").is_file())
            metrics = json.loads((output / "run_qc_metrics.json").read_text())
            self.assertEqual(metrics["overall"]["reads"], 3)
            self.assertEqual([row["barcode"] for row in metrics["barcodes"]], ["barcode01", "barcode02"])

    def test_combined_fastq_header_barcodes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fastq = Path(temporary) / "combined.fastq"
            self.write_fastq(
                fastq,
                [
                    ("read1 barcode=barcode01", "AAAA", "IIII"),
                    ("read2 BC:Z:barcode02", "CCCC", "IIII"),
                    ("read3", "GGGG", "IIII"),
                ],
            )
            stats = QC.analyse_fastqs(QC.discover_fastqs(fastq), "auto", 0, 0)
            self.assertEqual(stats["barcode01"].reads, 1)
            self.assertEqual(stats["barcode02"].reads, 1)
            self.assertEqual(stats["combined"].reads, 1)

    def test_malformed_fastq_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            fastq = Path(temporary) / "bad.fastq"
            fastq.write_text("@read1\nAAAA\n+\nIII\n", encoding="ascii")
            with self.assertRaisesRegex(ValueError, "length mismatch"):
                list(QC.read_fastq(fastq))


if __name__ == "__main__":
    unittest.main()
