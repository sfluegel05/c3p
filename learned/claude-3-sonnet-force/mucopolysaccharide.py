"""
Classifies: CHEBI:37395 mucopolysaccharide
"""
Based on the provided outcomes, the previous program appears to have some shortcomings in accurately classifying mucopolysaccharides. Here are my observations and suggestions for improvement:

1. **False Negatives**: The program missed several molecules that are known to be mucopolysaccharides, such as Pchaeglobosal B, Dendroamide B, and EH21A3. This suggests that the current patterns used for detecting uronic acid and glycosamine units may be too specific or restrictive.

2. **False Positives**: The program incorrectly classified Lanciferine as a mucopolysaccharide, even though it does not seem to contain the required structural features.

3. **Alternating Pattern Recognition**: While the program checks for the presence of uronic acid and glycosamine units, it does not explicitly verify if these units are arranged in an alternating pattern, which is a key feature of mucopolysaccharides.

4. **Esterification with Sulfuric Acid**: The program checks for the presence of sulfate groups, but it does not necessarily ensure that these sulfate groups are esterified to the polysaccharide backbone, as specified in the definition of mucopolysaccharides.

To improve the program's performance, I would suggest the following modifications:

1. **Refine Uronic Acid and Glycosamine Patterns**: Explore more general and flexible SMARTS patterns to capture a wider range of uronic acid and glycosamine units. You could consider using recursive SMARTS patterns or more advanced substructure matching techniques.

2. **Implement Alternating Pattern Check**: Develop a method to explicitly check if the uronic acid and glycosamine units are arranged in an alternating pattern along the polysaccharide backbone. This could involve analyzing the connectivity and sequence of the units.

3. **Ensure Sulfate Esterification**: In addition to checking for the presence of sulfate groups, verify that these sulfate groups are indeed esterified to the polysaccharide backbone. This could involve analyzing the connectivity and bonding patterns around the sulfate groups.

4. **Explore Additional Structural Features**: Investigate other structural features that are characteristic of mucopolysaccharides, such as the presence of specific glycosidic linkages, branching patterns, or specific substituents.

5. **Consider Machine Learning Approaches**: If the rule-based approach proves too challenging or inflexible, you could explore machine learning techniques, such as graph neural networks or fingerprint-based models, to learn the structural patterns of mucopolysaccharides from labeled data.

6. **Manually Review False Positives and Negatives**: Carefully review the false positives and false negatives, and use your chemical knowledge and intuition to determine if they are indeed misclassifications or if the benchmark data itself may contain errors or inconsistencies.

It's important to note that the classification of chemical entities can be a complex task, and achieving perfect accuracy may not always be possible, especially for broad and diverse chemical classes. However, by iteratively refining the patterns, incorporating additional structural features, and leveraging machine learning techniques if necessary, you can strive to improve the program's performance and robustness.