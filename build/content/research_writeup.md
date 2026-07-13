# Finding drug targets that act only when the immune system is switched on

## The problem

Many immune diseases, including rheumatoid arthritis, psoriasis, inflammatory
bowel disease, and multiple sclerosis, are driven by T cells that release too
much of particular signalling proteins called **cytokines**. Drugs that block
these cytokines already exist and work well, but they tend to blunt immune
signalling everywhere, all the time, which raises the risk of infection.

A more precise strategy is to target the immune response only *when it is
active*. A resting T cell and a stimulated one run very different genetic
programs. If we can find a gene that controls a disease-relevant cytokine
**only in the stimulated state, and not at rest**, then interfering with that
gene should calm an over-active response while leaving the resting immune
system largely untouched. Those genes, which we call context-specific
regulators, are what this project set out to find. It then asks a practical
follow-up question: of the genes we find, which ones could realistically be
turned into drugs, and by what kind of drug?

## The data

We built on a landmark public dataset from the Marson and Pritchard labs (Zhu,
Dann et al. 2025): a genome-scale **Perturb-seq** screen in primary human CD4+
T cells. In plain terms, the authors switched off essentially every gene, one
at a time, and read out the consequences on gene expression, in resting cells
and in cells stimulated for 8 and 48 hours. The full dataset is enormous
(hundreds of gigabytes). Rather than reprocess raw data, we used the authors'
precomputed, statistically rigorous effect sizes and streamed just the slice we
needed (the cytokine genes) directly from cloud storage, so the analysis runs
on a laptop.

## Finding context-specific regulators

The key statistical idea is to test an **interaction**, not two separate
comparisons. It is not enough for a gene to affect a cytokine in stimulated
cells; the effect has to be *specific* to that state. So for every gene, for
each of 30 cytokines, we tested whether its effect in stimulated cells is
genuinely different from its effect at rest, and corrected for the enormous
number of tests being run genome-wide.

A gene counts as a hit only if its knockdown changes a cytokine under
stimulation, has essentially no effect at rest, and passes that
interaction test. We then scored each gene for **selectivity**, rewarding
genes that focus on one or a few cytokines over "master switch" genes that
reshape the whole cell (which make messier drug targets), and applied a
reproducibility filter, keeping only genes whose knockdown was confirmed to hit
its intended target and to behave consistently across different human donors.

The result is a frozen shortlist of **128 high-confidence, context-specific
cytokine regulators**. Some are *suppressors* (switching them off lowers a
cytokine, giving anti-inflammatory leads); others are *inducers* (switching
them off raises a cytokine, giving leads for boosting immunity, for example in
cancer). Reassuringly,
known master regulators of T-cell identity land exactly where the selectivity
score predicts, which supports the method.

## From gene to drug: two ways to hit a target

A promising gene is only useful if it can actually be drugged. Historically,
the pharmaceutical industry has focused on proteins with a well-defined pocket
(for small-molecule pills) or on the cell surface (for antibodies), leaving a
large fraction of the genome branded "undruggable."

That framing is out of date. Genes that cannot be drugged at the protein level
are often perfectly targetable at the **RNA level**, using antisense
oligonucleotides (ASOs) or siRNA, nucleic-acid drugs that destroy a gene's
messenger RNA before it is ever made into protein. Several are already approved
medicines. So we scored all 128 genes on **two independent axes**: how druggable
each is as a protein (small molecule or antibody), and how druggable each is at
the RNA level (siRNA or ASO), the latter informed by transcript accessibility,
stability, and where in the cell the RNA resides.

Plotting both axes against each other turns "druggable vs undruggable" into a
map with a route for almost every gene. Of the 128, only a handful score poorly
on both axes; the great majority are addressable, and for each gene we
nominate a concrete modality. A meaningful group of otherwise hard-to-drug
regulators fall cleanly into the **RNA-preferred** category, where an
oligonucleotide strategy is the natural path forward.

## Why it matters

The output is a ranked, druggability-annotated shortlist of new candidate
targets for tuning specific cytokines in activated T cells, each paired with a
suggested therapeutic modality. It is a reproducible pipeline others can point
at a different cytokine, cell type, or disease context. And by scoring the RNA
axis as a first-class option rather than a fallback, it deliberately keeps in
play the many regulators that traditional target-discovery would have discarded
as undruggable, precisely the space where modern RNA therapeutics are
strongest.
