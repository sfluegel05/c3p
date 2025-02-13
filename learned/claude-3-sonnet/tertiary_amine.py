"""
Classifies: CHEBI:32876 tertiary amine
"""
The previous program attempts to classify a molecule as a tertiary amine by checking if it contains a nitrogen atom with exactly three carbon-containing substituents. While this approach is valid, it fails to account for some edge cases and additional requirements for tertiary amines, leading to the low F1 score.

Here are some potential issues and improvements:

1. **Missed cases with quaternary amines or N-oxides**: The program only checks for nitrogen atoms with three carbon substituents, but it should also exclude cases where the nitrogen is quaternary (four substituents) or an N-oxide.

2. **Handling aromatic rings**: The current implementation does not distinguish between aliphatic and aromatic carbon substituents. It may be necessary to treat aromatic rings as a single substituent rather than counting each carbon atom separately.

3. **Checking for hydrocarbyl groups**: The definition of a tertiary amine requires the substituents to be hydrocarbyl groups, which means they should only contain carbon and hydrogen atoms. The current implementation does not check for this condition.

4. **Handling charged and zwitterionic species**: Some of the false positives may be due to the program incorrectly identifying charged or zwitterionic species as tertiary amines.

5. **Considering structural context**: In some cases, the program may need to consider the broader structural context of the molecule to accurately identify tertiary amines. For example, it may need to account for specific functional groups or ring systems.

To improve the program, you could consider the following steps:

1. **Exclude quaternary amines and N-oxides**: Check the number of substituents on the nitrogen atom and exclude cases where there are more than three or where the nitrogen is an N-oxide.

2. **Handle aromatic rings**: Implement a method to identify aromatic rings and treat them as a single substituent when counting carbon substituents on the nitrogen atom.

3. **Check for hydrocarbyl groups**: Ensure that the substituents on the nitrogen atom only contain carbon and hydrogen atoms by checking the atomic numbers of the neighboring atoms.

4. **Handle charged and zwitterionic species**: Implement a method to identify charged or zwitterionic species and exclude them from being classified as tertiary amines.

5. **Consider structural context**: Incorporate additional checks or rules based on specific functional groups, ring systems, or other structural features that may help identify tertiary amines more accurately.

6. **Expand the training set**: Obtain a larger and more diverse set of positive and negative examples to improve the program's ability to generalize and handle edge cases.

By addressing these issues and incorporating the necessary improvements, you should be able to enhance the accuracy and robustness of the tertiary amine classification program.