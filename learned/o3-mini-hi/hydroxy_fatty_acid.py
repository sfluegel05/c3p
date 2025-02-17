"""
Classifies: CHEBI:24654 hydroxy fatty acid
"""
#!/usr/bin/env python
"""
Classifies: Hydroxy Fatty Acid
Definition:
  A molecule is considered a hydroxy fatty acid if it:
    - Contains a terminal carboxylic acid group of the form C(=O)[OX2H] 
      where the acid carbon is attached to exactly one other carbon.
    - Contains at least one extra –OH substituent that is not part of the acid group.
      (Extra –OH means an oxygen with a hydrogen that is bound to a carbon that is not
       a carbonyl-carbon.)
    - Has a minimum of 4 carbon atoms.
    - For molecules with 8 or more carbons:
         * The fraction of heavy atoms (non‐H) that are carbon is at least 0.75
         * At least 85% of the carbon atoms are non‐aromatic.
    - Contains at most one ring (small rings like epoxides are allowed).
Molecules not meeting these criteria are not hydroxy fatty acids.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a hydroxy fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as a hydroxy fatty acid, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens to help identify OH groups properly.
    mol = Chem.AddHs(mol)
    
    # Check that there is at most one ring (small epoxide rings allowed)
    if mol.GetRingInfo().NumRings() > 1:
        return False, "Molecule has too many rings to be a fatty acid"
    
    # --- Identify the terminal carboxylic acid group ---
    # Use a SMARTS pattern for a non-ring carboxylic acid group: [C;!R](=O)[O;H1]
    acid_pattern = Chem.MolFromSmarts("[C;!R](=O)[O;H1]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid group found (not a fatty acid)"
    
    valid_acid_match = None
    for match in acid_matches:
        # match[0] is the acid carbon, match[1] is the hydroxyl oxygen in the acid group.
        acid_carbon = mol.GetAtomWithIdx(match[0])
        # Terminal acid: the acid carbon must be attached to exactly one carbon.
        carbon_neighbors = [nbr for nbr in acid_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) == 1:
            valid_acid_match = match
            break
    if valid_acid_match is None:
        return False, "No terminal carboxylic acid group found (acid carbon not attached to exactly one carbon)"
    acid_oh_idx = valid_acid_match[1]  # to exclude the acid -OH later
    
    # --- Count carbon atoms and assess aliphaticity ---
    carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    n_carbons = len(carbons)
    if n_carbons < 4:
        return False, f"Too few carbons ({n_carbons}) to be a fatty acid"
    
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() != 1]
    n_heavy = len(heavy_atoms)
    if n_carbons >= 8 and n_heavy > 0:
        n_heavy_carbon = sum(1 for atom in heavy_atoms if atom.GetAtomicNum() == 6)
        if (n_heavy_carbon / n_heavy) < 0.75:
            return False, "Molecule is not sufficiently aliphatic (low fraction of carbon among heavy atoms)"
        n_aliphatic_carbons = sum(1 for atom in carbons if not atom.GetIsAromatic())
        if (n_aliphatic_carbons / n_carbons) < 0.85:
            return False, "Molecule is not sufficiently aliphatic (too many aromatic carbons)"
    
    # --- Identify extra hydroxyl (-OH) substituents outside the acid group ---
    extra_oh_found = False
    # Loop over all oxygen atoms in the molecule.
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "O":
            continue
        # Check that the oxygen has exactly one explicitly attached hydrogen.
        h_count = sum(1 for nbr in atom.GetNeighbors() if nbr.GetSymbol() == "H")
        if h_count != 1:
            continue
        # Exclude the hydroxyl oxygen that belongs to the carboxylic acid group.
        if atom.GetIdx() == acid_oh_idx:
            continue
        
        # Identify a carbon neighbor (an -OH group should be bound to a carbon)
        carbon_neighbor = None
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                carbon_neighbor = nbr
                break
        if carbon_neighbor is None:
            continue
        
        # Check that the carbon neighbor is not a carbonyl center.
        is_carbonyl = False
        for nbr in carbon_neighbor.GetNeighbors():
            if nbr.GetAtomicNum() == 8 and nbr.GetIdx() != atom.GetIdx():
                bond = mol.GetBondBetweenAtoms(carbon_neighbor.GetIdx(), nbr.GetIdx())
                if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
                    is_carbonyl = True
                    break
        if is_carbonyl:
            continue
        
        # Found an extra hydroxyl substituent.
        extra_oh_found = True
        break
    if not extra_oh_found:
        return False, "No hydroxy substituent found outside the carboxylic acid group"
    
    return True, ("Molecule contains a terminal carboxylic acid group, " +
                  ("a long aliphatic chain, " if n_carbons >= 8 else "") +
                  "and extra hydroxy substituent(s), classifying it as a hydroxy fatty acid")

# Example usage:
if __name__ == '__main__':
    test_smiles = [
        # True positives (should be classified as hydroxy fatty acids)
        "C(C(O)=O)C[C@@H](/C=C/C=C\\C/C=C\\C/C=C\\C=C\\[C@H](C/C=C\\CC)O)O",  # resolvin D6
        "OC(=CC)C(O)=O",  # 2-hydroxy-2-butenoic acid
        "C(C(O)=O)C/C=C\\C[C@H](\\C=C\\C=C/C=C/C=C/[C@H]([C@@H](C/C=C\\CC)O)O)O",  # aspirin-triggered resolvin D2
        "CCCC\\C=C/C[C@@H]1O[C@H]1\\C=C\\[C@H](O)C\\C=C/CCCC(O)=O",  # (8R)-hydroxy-(11S,12S)-epoxyicosa-(5Z,9E,14Z)-trienoic acid
        "CC[C@@H](O)CCCCCCCCCCCCCCCC[C@@H](O)CC(O)=O",  # (3R,20R)-3,20-dihydroxyhenicosanoic acid
        "O[C@@H]([C@H](C)C(O)=O)C",  # 3-Hydroxy-2-methyl-[R-(R,S)]-butanoic acid
        # False positive example (should be rejected)
        "OC(=O)CCC(CCCC)CC",  # 4-Ethyloctanoic acid: lacks any extra –OH
    ]
    
    for s in test_smiles:
        res, reason = is_hydroxy_fatty_acid(s)
        print(f"SMILES: {s}\nResult: {res}, Reason: {reason}\n")