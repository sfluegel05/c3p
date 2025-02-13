"""
Classifies: CHEBI:35267 quaternary ammonium ion
"""
The previous program attempted to classify quaternary ammonium ions by checking for the presence of a quaternary nitrogen atom with four univalent substituents and additional checks on molecular weight, number of atoms, and number of bonds. However, the approach has some limitations:

1. **False Positives**: The program incorrectly classifies molecules containing a quaternary nitrogen atom but not belonging to the quaternary ammonium ion class. Examples include carnitines, indole derivatives, and other compounds containing a quaternary nitrogen but not specifically quaternary ammonium ions.

2. **False Negatives**: The program misses some valid quaternary ammonium ions, such as simple quaternary ammonium salts like pentamethonium and choline derivatives without a phosphocholine group.

3. **Overly Restrictive Filters**: The additional checks on molecular weight, number of atoms, and number of bonds may be too restrictive and exclude valid quaternary ammonium ions that fall outside the specified ranges.

To improve the program, we can consider the following strategies:

1. **Use a More Specific SMARTS Pattern**: Instead of checking for any quaternary nitrogen atom, we can use a more specific SMARTS pattern that captures the core structure of quaternary ammonium ions. This can be achieved by incorporating the nitrogen atom and the four substituents in the pattern.

2. **Handle Simple Quaternary Ammonium Salts**: Include a separate check for simple quaternary ammonium salts, such as pentamethonium and choline derivatives, which may not match the more complex SMARTS pattern.

3. **Remove Overly Restrictive Filters**: Consider removing or relaxing the additional filters on molecular weight, number of atoms, and number of bonds, as they may be excluding valid quaternary ammonium ions.

4. **Incorporate Negative Filtering**: Use additional SMARTS patterns or chemical rules to filter out known false positives, such as carnitines, indole derivatives, and other common compounds that may be misclassified.

5. **Consult Chemical Databases or Experts**: Consult chemical databases or domain experts to identify more comprehensive rules or patterns for accurately classifying quaternary ammonium ions, as the class can have diverse structures and substituents.

By incorporating these strategies, the program can potentially improve its accuracy in classifying quaternary ammonium ions while reducing false positives and false negatives.