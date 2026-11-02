import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks/public_modeling"))
from paired_predictions import compare


class PairedPredictionsTests(unittest.TestCase):
    def setUp(self):
        self.rows = [
            dict(model=model, donor=str(i), fold="holdout", observed=i, predicted=i+error)
            for model, error in [("RNA", 2), ("combined", 1)] for i in range(3)
        ]

    def test_constant_donor_gain(self):
        result = compare(self.rows, "RNA")
        self.assertEqual(result[1], ("combined", 3, 1, 1, 1, 1))

    def test_incomplete_or_repeated_donors(self):
        for rows in (self.rows[:-1], self.rows + self.rows[:1]):
            with self.assertRaises(ValueError):
                compare(rows, "RNA")

    def test_fold_outcome_and_finite_predictions(self):
        for update in (dict(fold="training"), dict(observed=4), dict(predicted=float("nan"))):
            rows = self.rows[:-1] + [{**self.rows[-1], **update}]
            with self.assertRaises(ValueError):
                compare(rows, "RNA")
