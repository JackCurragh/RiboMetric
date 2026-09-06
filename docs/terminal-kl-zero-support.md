# Terminal KL divergence when background support is zero

The raw terminal dinucleotide KL function returns positive infinity if a
positive observed probability has zero expected probability. Zero observed
mass contributes zero, including a 0/0 term. Missing expected keys still raise
an input error. Negative, non-finite and greater-than-one probabilities raise
ValueError rather than silently enter the score.

The existing goodness transformation `1/(1 + raw_KL)` maps infinity to zero.
For example, observed AA=1, CC=0 and expected AA=0, CC=1 must not receive a
perfect-agreement score.

The exported metrics keep the existing score keys. For raw exports:

- finite divergence: the raw key is a number and its `_raw_status` is `finite`;
- zero-background support mismatch: the raw key is null and its `_raw_status`
  is `zero_background_support`.

Both 5prime and 3prime statuses are explicit. A JSON null in this case means
infinite raw divergence, not a zero divergence or a failed/omitted calculation.
This avoids the non-standard JSON Infinity token while preserving the reason.
The Python raw-divergence function itself continues to return a float.

This patch does not invent a smoothing policy. An empirical zero is not proof
of biological impossibility; a future background-estimation/smoothing change
needs its own scientific rationale. No claim is made about the frequency of
this input in BAM/FASTA datasets or effects on published results.

Reference for the zero-support convention:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.rel_entr.html
