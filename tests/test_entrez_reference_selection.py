import unittest
from unittest.mock import patch, MagicMock
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os, sys

# Ensure we import the local package from src/ without installation
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src'))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from wildtype_aligner.engine import AlignmentEngine


class TestEntrezReferenceSelection(unittest.TestCase):
    @patch('wildtype_aligner.engine.SeqIO.read')
    @patch('wildtype_aligner.engine.Entrez.efetch')
    @patch('wildtype_aligner.engine.Entrez.read')
    @patch('wildtype_aligner.engine.Entrez.esummary')
    @patch('wildtype_aligner.engine.Entrez.esearch')
    def test_prefers_refseq_np_accession(self, mock_esearch, mock_esummary, mock_read, mock_efetch, mock_seqio_read):
        # esearch returns IdList with several IDs
        mock_esearch.return_value = MagicMock(close=MagicMock())
        search_record = {"IdList": ["111", "222", "333"]}

        # esummary returns docs where one has NP_ accession (choose that ID)
        mock_esummary.return_value = MagicMock(close=MagicMock())
        summary_record = {
            "DocumentSummarySet": {
                "DocumentSummary": [
                    {"Id": "111", "AccessionVersion": "XP_0001"},
                    {"Id": "222", "AccessionVersion": "NP_012345"},
                    {"Id": "333", "AccessionVersion": "YP_999999"},
                ]
            }
        }

        # read is called twice: once for esearch, once for esummary
        mock_read.side_effect = [search_record, summary_record]

        # efetch + SeqIO.read returns a sequence
        mock_efetch.return_value = MagicMock(close=MagicMock())
        ref_seq = SeqRecord(Seq("MKTIIALSYIFCLVFAD"), id="NP_012345")
        mock_seqio_read.return_value = ref_seq

        engine = AlignmentEngine(max_retries=2)
        seqrec = engine._fetch_sepi_reference('acrB', 'Escherichia', 'coli')

        self.assertIn('NP_', seqrec.id)
        self.assertEqual(str(seqrec.seq), "MKTIIALSYIFCLVFAD")


if __name__ == '__main__':
    unittest.main()
