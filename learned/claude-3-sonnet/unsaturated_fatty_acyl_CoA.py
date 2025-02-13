"""
Classifies: CHEBI:51006 unsaturated fatty acyl-CoA
"""
The previous program had a few issues that led to the poor performance:

1. **Overly Strict Pattern Matching**: The program relied heavily on SMARTS pattern matching, which can be too strict and miss valid structures. For example, the CoA backbone pattern was very specific and would not match slightly different representations of the same backbone.

2. **Inappropriate Filters**: Some of the filters were too strict or not appropriate for the class. For example, requiring a minimum of 8 rotatable bonds and 12 carbon atoms would exclude smaller unsaturated fatty acyl-CoAs, which are still valid members of the class.

3. **Molecular Weight Filter**: The molecular weight filter (> 500 Da) was too high and would exclude many valid unsaturated fatty acyl-CoAs, which typically have lower molecular weights.

4. **Lack of Stereochemistry Considerations**: The program did not account for stereochemistry, which is essential for correctly identifying unsaturated fatty acyl-CoAs with different double bond configurations.

To improve the program, we should:

1. **Use More Flexible Matching**: Instead of relying solely on strict SMARTS patterns, we should use more flexible methods to identify key functional groups and substructures.

2. **Adjust Filters**: Relax or remove filters that are too strict or inappropriate for the class, such as the minimum rotatable bond and carbon count filters.

3. **Adjust Molecular Weight Range**: Use a more appropriate molecular weight range for unsaturated fatty acyl-CoAs, which is typically lower than 500 Da.

4. **Consider Stereochemistry**: Incorporate stereochemistry checks to ensure that the double bond configurations match those expected for unsaturated fatty acyl-CoAs.

5. **Use Additional Structural Features**: Exploit other structural features characteristic of unsaturated fatty acyl-CoAs, such as the presence of a thioester bond between the fatty acid and CoA moieties, and the presence of an unsaturated carbon chain with a specific number of double bonds.

By addressing these issues, we can develop a more robust and accurate program for classifying unsaturated fatty acyl-CoAs based on their SMILES strings.