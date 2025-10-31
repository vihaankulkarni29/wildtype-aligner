import unittest
from unittest.mock import patch
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os, sys

# Ensure we import the local package from src/ without installation
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src'))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from wildtype_aligner.engine import AlignmentEngine


class TestBatchAlignmentMocked(unittest.TestCase):
    def make_queries(self):
        # IDs encode gene as last underscore-part
        q1 = SeqRecord(Seq("MKTIIALSYIFCLVFAD"), id="ACC1_acrA", description="")
        q2 = SeqRecord(Seq("MKTIIALSYIFCLVFAD"), id="ACC2_acrA", description="")
        q3 = SeqRecord(Seq("MKKLLTAAALVALATASA"), id="ACC3_acrB", description="")
        return [q1, q2, q3]

    @patch.object(AlignmentEngine, '_get_reference_with_cache')
    def test_batch_alignments_with_mock_refs(self, mock_get_ref):
        # Provide reference sequences per gene
        ref_map = {
            'acra': SeqRecord(Seq("MKTIIALSYIFCLVFAD"), id="REF_acrA"),
            'acrb': SeqRecord(Seq("MKKLLTAAALVALATASA"), id="REF_acrB"),
        }

        def ref_side_effect(gene, genus, species):
            return ref_map[gene.lower()]

        mock_get_ref.side_effect = ref_side_effect

        engine = AlignmentEngine()
        queries = self.make_queries()
        report, summary = engine.run_batch_alignments(queries, 'Escherichia', 'coli', show_progress=False)

        self.assertIsInstance(report, str)
        self.assertEqual(summary.total, 3)
        self.assertEqual(summary.failed, 0)
        self.assertEqual(summary.successful, 3)
        self.assertIn('# 1', report)
        self.assertIn('REF_acrA', report)
        self.assertIn('REF_acrB', report)

    @patch.object(AlignmentEngine, '_get_reference_with_cache')
    @patch.object(AlignmentEngine, '_align_single')
    def test_error_handling_counts_failures(self, mock_align, mock_get_ref):
        # Always return a valid ref
        ref = SeqRecord(Seq("MKTIIALSYIFCLVFAD"), id="REF")
        mock_get_ref.return_value = ref

        # First alignment raises, second returns None, third succeeds
        class DummyAlign:
            def __init__(self):
                self.score = 1
            def __getitem__(self, i):
                return "MKTIIALSYIFCLVFAD"

        mock_align.side_effect = [RuntimeError("boom"), None, DummyAlign()]

        engine = AlignmentEngine()
        queries = [
            SeqRecord(Seq("A"), id="Q1_acrA"),
            SeqRecord(Seq("A"), id="Q2_acrA"),
            SeqRecord(Seq("A"), id="Q3_acrA"),
        ]

        report, summary = engine.run_batch_alignments(queries, 'Escherichia', 'coli', show_progress=False)

        self.assertEqual(summary.total, 3)
        self.assertEqual(summary.successful, 1)
        self.assertEqual(summary.failed, 2)
        self.assertIn('Error aligning', report)
        self.assertIn('No alignment found', report)


if __name__ == '__main__':
    unittest.main()
