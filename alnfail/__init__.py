"""alnfail: an independent benchmark of read aligners with known ground truth.

Project 26, "Where Do Aligners Actually Fail?" (AITU, Introduction to Bioinformatics).

The package is organised by task:

    coords      coordinate conventions and contig-name mapping        (task 2)
    fetch       scripted, logged retrieval of references and reads    (task 2)
    annot       repeat / duplication / feature annotations -> strata  (tasks 2, 7)
    variants    truth-set variants applied to simulated reads         (tasks 2, 7)
    qc          read QC with quantified filtering decisions           (NGS / TGS)
    profile     learn a platform error profile from real reads        (task 3)
    simulate    simulate reads whose true origin is known             (tasks 4, 5)
    samtable    alignments -> columnar table (Parquet)                (task 6)
    evaluate    correct / misplaced / unmapped per read               (tasks 4-6)
    summarise   accuracy, MAPQ calibration, stratified failures       (tasks 6, 7)
    consensus   real-data validation without per-read truth           (task 8)
    mapq_theory what MAPQ is supposed to mean                         (task 1)
"""

__version__ = "1.0.0"
