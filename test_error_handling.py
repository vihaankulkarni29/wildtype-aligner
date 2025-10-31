"""
Test script to verify error handling functionality
Creates intentionally problematic sequences to test error tracking
"""
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from src.wildtype_aligner.engine import AlignmentEngine

# Create test sequences with intentional issues
test_sequences = [
    # Valid sequences
    SeqRecord(Seq("MKLPQRSTV"), id="valid_seq_1", description="gene=acrA protein test"),
    SeqRecord(Seq("MLKRSTPQV"), id="valid_seq_2", description="gene=acrA protein test"),
    
    # Empty sequence (should trigger error)
    SeqRecord(Seq(""), id="empty_seq", description="gene=acrA empty sequence"),
    
    # Very short sequence (might cause alignment issues)
    SeqRecord(Seq("MK"), id="short_seq", description="gene=acrA short protein"),
    
    # Another valid one
    SeqRecord(Seq("MKLPQRSTVAA"), id="valid_seq_3", description="gene=acrA protein test"),
]

print("="*70)
print("Testing Error Handling with Mixed Valid/Invalid Sequences")
print("="*70)

engine = AlignmentEngine()

print("\n🧪 Running batch alignments with error tracking...")
try:
    report, summary = engine.run_batch_alignments(
        test_sequences,
        genus="Escherichia",
        species="coli",
        error_log="test_error.log",
        show_progress=True
    )
    
    print("\n" + "="*70)
    print("RESULTS:")
    print("="*70)
    print(summary.get_summary())
    
    if summary.failed > 0:
        print(f"\n✅ Error log saved to: test_error.log")
        print("\nError details:")
        for seq_id, error_msg in summary.errors:
            print(f"  - {seq_id}: {error_msg}")
    else:
        print("\n✅ All sequences aligned successfully!")
    
    print("\n" + "="*70)
    print("TEST COMPLETE")
    print("="*70)
    
except Exception as e:
    print(f"\n❌ Test failed with error: {e}")
    import traceback
    traceback.print_exc()
