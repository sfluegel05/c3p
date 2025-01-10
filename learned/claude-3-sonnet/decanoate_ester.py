"""
Classifies: CHEBI:87658 decanoate ester
"""
Looking at the outcomes, the main issue is that the program has too many false positives - it's identifying many compounds as decanoate esters that aren't actually decanoate esters. The key problems appear to be:

1. The SMARTS pattern for the decanoyl chain (`[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]`) is too flexible with the `~` bonds, allowing for double bonds and rings within the chain.

2. The chain validation logic isn't strict enough - it allows for branched chains and substitutions that shouldn't be part of a decanoate ester.

3. We need to ensure the chain is specifically 10 carbons long, not longer.

Here's how we can improve:

1. Use more specific SMARTS patterns:
   - Use `-` instead of `~` to specify single bonds
   - Add explicit hydrogens where needed
   - Specify saturated carbons (CX4) for most of the chain

2. Add stricter chain validation:
   - Count exact number of carbons
   - Ensure linear chain structure
   - Verify no branching except at the end

3. Add additional checks:
   - Verify the carbon chain is saturated (except for the carbonyl carbon)
   - Check that the ester oxygen is not part of a larger group like phosphate

Here's the reasoning for a new implementation:
1. First check for basic ester group
2. Find all potential decanoyl chains
3. For each chain:
   - Verify it's exactly 10 carbons
   - Verify it's linear (no branches except at ends)
   - Verify it's saturated (all single bonds)
   - Verify proper connectivity to ester group
4. Additional validation to exclude false positives:
   - Check bond types
   - Verify chain length
   - Ensure proper connectivity

This should help reduce false positives while maintaining true positives.

Would you like me to implement this improved version?