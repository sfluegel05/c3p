"""
Classifies: CHEBI:17962 CDP-diacylglycerol
"""
Based on the failed examples and the provided definition of CDP-diacylglycerol, it seems that the previous program had a few issues:

1. **CDP-glycerol backbone pattern**: The SMARTS pattern used for the CDP-glycerol backbone was too specific and did not account for the variations in the cytidine diphosphate moiety. The pattern should be more general to cover different tautomeric forms and protonation states.

2. **Ester group count**: The program assumed exactly two ester groups, which is not always the case. Some examples had additional ester or hydroxyl groups, leading to false negatives.

3. **Fatty acid chain detection**: The SMARTS pattern used for detecting fatty acid chains was too simplistic and did not account for various substituents or unsaturations present in the examples.

4. **Oxygen count**: The program assumed a fixed number of oxygens (15), which is incorrect as the number can vary based on the acyl chains and additional functional groups.

To improve the program, we can make the following changes:

1. **CDP-glycerol backbone pattern**: Use a more general SMARTS pattern for the CDP-glycerol backbone, such as `"[C@H]1[C@@H]([C@H]([C@@H](O1)N1C=NC(=NC1=O)N)O)OP(=O)(O)OP(=O)(O)O"`. This pattern matches the core glycerol moiety with the cytidine diphosphate group, without being too specific about the protonation states or tautomeric forms.

2. **Ester group count**: Instead of checking for a fixed number of ester groups, we can check if there are at least two ester groups present, allowing for additional functional groups.

3. **Fatty acid chain detection**: Use a more flexible SMARTS pattern that accounts for unsaturations, substituents, and various chain lengths, such as `"[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]"`.

4. **Oxygen count**: Remove the check for a specific oxygen count, as it can vary based on the acyl chains and additional functional groups.

5. **Additional checks**: Optionally, we can add checks for the presence of long aliphatic chains (e.g., molecular weight, rotatable bond count) to further refine the classification.

By incorporating these changes, the program should be able to correctly classify the provided examples as CDP-diacylglycerols, while also being more robust to variations in the structure.