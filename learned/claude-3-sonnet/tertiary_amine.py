"""
Classifies: CHEBI:32876 tertiary amine
"""
Based on the outcomes and the previous code, it seems that the program is not accurately identifying tertiary amines in all cases. Here are some potential issues and ways to improve the program:

1. **Handling aromatic systems**: The current code checks for aromatic substituents on the nitrogen atom, but it doesn't handle aromatic systems correctly. For example, in `triflupromazine` (SMILES: `CN(C)CCCN1c2ccccc2Sc2ccc(cc12)C(F)(F)F`), the program fails to recognize the aromatic system as a hydrocarbyl substituent. To address this, the code could use a more robust method to identify aromatic systems, such as the `GetIsAromatic` function from the `rdkit.Chem.rdchem` module.

2. **Considering bond orders**: The current code does not consider bond orders when determining if a substituent is a hydrocarbyl group. This could lead to incorrect classifications for molecules with double or triple bonds in the substituents. For example, in `Diallat` (SMILES: `CC(C)N(C(C)C)C(=O)SC\C(Cl)=C\Cl`), the program fails to recognize the double bond in the substituent. To improve this, the code could check the bond orders between atoms in the substituent.

3. **Handling cyclic structures**: The current code does not handle cyclic structures correctly. For example, in `dicyclomine` (SMILES: `CCN(CC)CCOC(=O)C1(CCCCC1)C1CCCCC1`), the program fails to recognize the cyclic substituents as hydrocarbyl groups. To address this, the code could use a recursive function to traverse the molecule and identify cyclic structures.

4. **Checking for specific functional groups**: The current code does not check for specific functional groups that should be excluded from the definition of a tertiary amine. For example, in `(6S)-6-hydroxyhyoscyamine` (SMILES: `CN1[C@H]2C[C@@H](C[C@@H]1[C@@H](O)C2)OC(=O)C(CO)c1ccccc1`), the program incorrectly classifies the molecule as a tertiary amine due to the presence of an alcohol group. To address this, the code could check for specific functional groups and exclude molecules containing them.

5. **Handling charged species**: The current code does not handle charged species correctly. For example, in `malachite green isothiocyanate cation` (SMILES: `CN(C)c1ccc(cc1)C(c1ccc(cc1)N=C=S)=C1C=CC(C=C1)=[N+](C)C`), the program fails to recognize the quaternary nitrogen atom. To address this, the code could check for charged atoms and handle them appropriately.

To improve the program, you could implement the following changes:

1. Use a more robust method to identify aromatic systems, such as the `GetIsAromatic` function from the `rdkit.Chem.rdchem` module.
2. Check bond orders between atoms in the substituents to ensure they are hydrocarbyl groups.
3. Implement a recursive function to traverse the molecule and identify cyclic structures.
4. Check for specific functional groups that should be excluded from the definition of a tertiary amine.
5. Handle charged species correctly by checking for charged atoms and treating them appropriately.

By addressing these issues, the program should be able to more accurately classify tertiary amines based on their SMILES strings.