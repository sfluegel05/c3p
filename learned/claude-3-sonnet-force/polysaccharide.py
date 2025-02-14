"""
Classifies: CHEBI:18154 polysaccharide
"""
Analyzing the previous program and the outcomes, it seems that the main issue is that the program is too restrictive in its definition of a polysaccharide. The program checks for the presence of carbohydrate rings, glycosidic bonds, and multiple monosaccharide units, but some of the examples provided by the benchmark are more complex polysaccharides with additional modifications and substituents.

Here are some observations and potential improvements:

1. **Handling Complex Structures**: The program assumes a relatively simple structure of repeating monosaccharide units linked by glycosidic bonds. However, some of the false negatives seem to be more complex polysaccharides with additional substituents or modifications. To handle these cases, the program could be modified to allow for a broader range of substituents and modifications while still checking for the core polysaccharide structure.

2. **Relaxing Molecular Weight Constraint**: The molecular weight constraint of >500 Da is too restrictive for some polysaccharides. While it is generally true that polysaccharides have a high molecular weight, there can be exceptions, especially for smaller polysaccharides or those with substituents that lower the molecular weight. Relaxing or removing this constraint could improve the program's performance.

3. **Considering Additional Patterns**: The program currently checks for specific patterns of carbohydrate rings and glycosidic bonds. However, there may be other patterns or substructures that are characteristic of polysaccharides that could be incorporated into the program. Analyzing the false negatives and false positives could provide insights into additional patterns to consider.

4. **Handling Outliers**: As mentioned, there may be occasional and systematic mistakes in the benchmark used for evaluation. In such cases, if the program's classifications align with your understanding of the chemical class, it may be reasonable to ignore outliers and focus on improving the program's performance for the majority of cases. However, this should be done with caution and proper reasoning.

5. **Incorporating Domain Knowledge**: While the program relies on structural patterns, incorporating domain knowledge about polysaccharides could help improve the classifications. For example, considering the types of monosaccharide units commonly found in polysaccharides, the typical linkage patterns, or other biochemical properties could provide additional criteria for classification.

Overall, the program's current approach is a good starting point, but it may need to be more flexible and incorporate additional patterns and domain knowledge to handle the diverse range of polysaccharide structures encountered in the benchmark. Iterative refinement based on analyzing false positives and false negatives would be a reasonable approach to improve the program's performance.