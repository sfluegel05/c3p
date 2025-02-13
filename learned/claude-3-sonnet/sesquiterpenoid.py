"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
The previous program had a few issues that led to the low F1 score:

1. **Incomplete sesquiterpene skeleton patterns**: The list of SMARTS patterns for identifying sesquiterpene skeletons was not comprehensive enough. While it covered some common skeletons, it missed many others, leading to false negatives.

2. **Over-reliance on skeleton patterns**: The program heavily relied on finding a specific skeleton pattern. If none were found, it then checked other criteria like ring count and rotatable bonds. However, these criteria alone are not sufficient to definitively classify a molecule as a sesquiterpenoid, leading to false positives.

3. **Strict carbon count requirement**: The program required exactly 15 carbon atoms for a sesquiterpenoid. While this is generally true, there are exceptions where sesquiterpenoids can have fewer or more carbon atoms due to rearrangements or modifications.

To improve the program, we can consider the following strategies:

1. **Expand the list of sesquiterpene skeleton patterns**: Manually curate a more comprehensive list of SMARTS patterns that cover a wider range of sesquiterpene skeletons. This will help reduce false negatives.

2. **Incorporate more specific structural features**: In addition to skeleton patterns, look for other structural features that are characteristic of sesquiterpenoids, such as the presence of specific functional groups (e.g., alcohols, esters, ethers), specific ring systems (e.g., decalins, bicyclo[3.3.1] systems), or specific carbon skeletal rearrangements.

3. **Relax the strict carbon count requirement**: Instead of requiring exactly 15 carbon atoms, allow for a range of carbon counts (e.g., 13-17) to account for potential modifications or rearrangements.

4. **Consider using machine learning models**: If curating a comprehensive set of rules and patterns becomes too complex, consider training a machine learning model (e.g., random forest, neural network) on a large dataset of labeled sesquiterpenoids and non-sesquiterpenoids. This could potentially learn the complex patterns more effectively than hand-crafted rules.

5. **Incorporate additional descriptors**: Explore the use of additional molecular descriptors (e.g., topological indices, molecular fingerprints) that could help differentiate sesquiterpenoids from other molecules.

By addressing these issues, the program's performance in classifying sesquiterpenoids should improve significantly.