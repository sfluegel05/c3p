"""
Classifies: CHEBI:13248 anilide
"""
Let me analyze the previous program and the outcomes:

Key issues identified:

1. False Negatives (Missed Anilides):
- The main issue is that the program missed several fentanyl derivatives and similar compounds where the nitrogen is tertiary (N-substituted) rather than secondary
- The current SMARTS pattern '[c][NH][CX3](=[OX1])' only matches secondary amides (with NH)
- We need to also match tertiary amides where the nitrogen is substituted

2. False Positives:
- Many of the false positives involve urea groups (-NH-C(=O)-NH-) being misidentified as amides
- Some involve more complex ring systems where the aromatic ring is part of a larger fused system
- Several cases involve carbamates (-NH-C(=O)-O-) being misidentified

Here's how we can improve the code:

1. Modify the SMARTS pattern to:
   - Match both secondary and tertiary N-phenyl amides
   - Exclude urea groups by ensuring the carbon isn't bonded to another nitrogen
   - Exclude carbamates by ensuring the carbon isn't bonded to oxygen

2. Add additional validation:
   - Better verification of the phenyl ring being a simple benzene ring
   - Check that the amide carbon isn't part of a urea or carbamate group

3. Improve ring detection:
   - Add checks to ensure we're matching a proper phenyl ring
   - Exclude cases where the ring is part of a larger fused system

Let me implement these improvements in the next code block.