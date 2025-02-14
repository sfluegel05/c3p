"""
Classifies: CHEBI:35346 11beta-hydroxy steroid
"""
Based on the provided outcomes, it appears that the previous program failed to correctly classify any of the positive or negative examples for the 11beta-hydroxy steroid class. This could be due to several reasons:

1. **Steroid backbone pattern**: The SMARTS pattern used to identify the steroid backbone might be too specific or restrictive. The pattern assumes a particular arrangement of rings and substituents, which may not cover all possible variations within the class. Additionally, the pattern includes an explicit hydrogen count, which can be problematic in SMILES representations.

2. **11beta-hydroxy group identification**: The SMARTS pattern used to identify the 11beta-hydroxy group might be too restrictive or not capturing all possible conformations and environments. The pattern assumes a specific ring context, which may not hold true for all molecules in the class.

3. **Additional features**: The additional checks for the number of ketones, double bonds, and molecular weight range might be too stringent or not representative of the entire class. These features can vary among different subclasses or derivatives of 11beta-hydroxy steroids.

To improve the program, we can consider the following approaches:

1. **Refine the steroid backbone pattern**: Analyze the positive examples and identify common structural patterns for the steroid backbone. Develop a more flexible SMARTS pattern that captures the essential ring system and connectivity while allowing for variations in substituents and stereochemistry.

2. **Improve the 11beta-hydroxy group identification**: Analyze the positive examples and identify common patterns or environments for the 11beta-hydroxy group. Develop a more flexible SMARTS pattern that can capture different conformations and ring contexts.

3. **Relax or adjust additional feature checks**: Analyze the positive and negative examples to understand the range of values for the number of ketones, double bonds, and molecular weights. Adjust the checks accordingly or consider removing them if they are not essential for the classification.

4. **Incorporate additional structural features**: Analyze the positive examples and identify any additional structural features or patterns that are characteristic of the 11beta-hydroxy steroid class. Incorporate these features into the classification logic using SMARTS patterns or other RDKit functionality.

5. **Use machine learning techniques**: If the class is challenging to define with explicit rules and patterns, consider using machine learning techniques, such as trained classifiers or fingerprint-based methods, to learn the structural features that distinguish the 11beta-hydroxy steroid class from other classes.

It's important to note that the chemical class definition provided ("Any 11-hydroxy steroid in which the hydroxy group at position 11 has beta- configuration.") might be incomplete or not accurately represent the entire class. In such cases, it may be necessary to consult additional sources or domain experts to refine the class definition and ensure that the program accurately captures the intended scope.

Additionally, if you have high confidence in some of the positive or negative examples based on your understanding of the chemical class, you could consider adjusting the benchmark or ignoring potential outliers while explaining your reasoning.