"""
Classifies: CHEBI:24654 hydroxy fatty acid
"""
#!/usr/bin/env python
"""
Classifies: Hydroxy Fatty Acid
Definition:
  A molecule is considered a hydroxy fatty acid if it:
    - Contains a terminal carboxylic acid group of the form C(=O)[O;H] 
      where the acid carbon is attached to exactly one other carbon.
    - Contains at least one extra –OH substituent (an O with one H attached)
      that is not the –OH of the acid group and is attached to a carbon that is not carbonyl.
    - Has a minimum of 4 carbon atoms.
    - If the molecule has 8 or more carbons then:
         * The fraction of heavy atoms (non‐H) that are carbon is high.
           (We require 0.70 for borderline cases (8–10 carbons) and 0.75 otherwise.)
         * At least 85% of the carbon atoms are non‐aromatic.
         * There is an unbroken chain of at least 4 non‐aromatic, non‐ring carbons.
    - If the molecule is short (<8 carbons) then it must be non‐aromatic.
    - Contains at most one ring (small rings such as epoxides are allowed).
Molecules failing these criteria are not hydroxy fatty acids.
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
    # Parse molecule and add explicit H (helps in functional group detection)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # Check that the molecule does not have too many rings.
    if mol.GetRingInfo().NumRings() > 1:
        return False, "Molecule has too many rings to be a fatty acid"
    
    # --- Identify the terminal carboxylic acid group ---
    # SMARTS: a non‐ring C(=O)[O;H] group. (The acid –OH is then explicitly identified.)
    acid_pattern = Chem.MolFromSmarts("[C;!R](=O)[O;H]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid group found (not a fatty acid)"
    
    valid_acid_match = None
    acid_oh_idx = None
    for match in acid_matches:
        acid_carbon = mol.GetAtomWithIdx(match[0])
        # Terminal acid: the acid carbon should be linked to exactly one carbon.
        carbon_neighbors = [nbr for nbr in acid_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) == 1:
            valid_acid_match = match
            acid_oh_idx = match[1]  # the index of the –OH oxygen in the acid group.
            break
    if valid_acid_match is None:
        return False, "No terminal carboxylic acid group found (acid carbon not attached to exactly one carbon)"
    
    # --- Basic atom count and aliphatic checks ---
    carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    n_carbons = len(carbons)
    if n_carbons < 4:
        return False, f"Too few carbons ({n_carbons}) to be a fatty acid"
    
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() != 1]
    n_heavy = len(heavy_atoms)
    
    # For molecules with fewer than 8 carbons, require that there is no aromatic carbon (to avoid aromatic acids)
    if n_carbons < 8:
        if any(atom.GetIsAromatic() for atom in carbons):
            return False, "Molecule is aromatic and too short to be a fatty acid"
    else:
        # For molecules with 8 or more carbons, require a high fraction of carbons.
        carbon_ratio = n_carbons / n_heavy if n_heavy > 0 else 0
        threshold = 0.70 if n_carbons <= 10 else 0.75
        if carbon_ratio < threshold:
            return False, f"Molecule is not sufficiently aliphatic (carbon/heavy-atom ratio {carbon_ratio:.2f} < {threshold})"
        # Also, at least 85% non-aromatic carbons should be present.
        non_aromatics = sum(1 for atom in carbons if not atom.GetIsAromatic())
        if (non_arom / n_carbons) < 0.85:
            return False, "Molecule is not sufficiently aliphatic (too many aromatic carbons)"
        # Require a long aliphatic chain: four or more consecutive non‐aromatic, non‐ring carbons.
        chain_pattern = Chem.MolFromSmarts("[C;!R;!$(c)]-[C;!R;!$(c)]-[C;!R;!$(c)]-[C;!R;!$(c)]")
        if not mol.HasSubstructMatch(chain_pattern):
            return False, "No long aliphatic chain found in a molecule with many carbons"
    
    # --- Identify extra hydroxyl (-OH) substituents outside the acid group ---
    extra_oh_found = False
    # Look for OH groups by iterating over oxygen atoms.
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "O":
            continue
        # Check that the oxygen has exactly one attached hydrogen.
        h_count = sum(1 for nbr in atom.GetNeighbors() if nbr.GetSymbol() == "H")
        if h_count != 1:
            continue
        # Skip the acid OH that belongs to the carboxylic acid group.
        if atom.GetIdx() == acid_oh_idx:
            continue
        # Verify that the oxygen is attached to a carbon.
        neighbor_carbons = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if not neighbor_carbons:
            continue
        # For each carbon neighbor, check it is not a carbonyl carbon.
        is_carbonyl_neighbor = False
        for carbon in neighbor_carbons:
            for nbr in carbon.GetNeighbors():
                if nbr.GetAtomicNum() == 8 and nbr.GetIdx() != atom.GetIdx():
                    bond = mol.GetBondBetweenAtoms(carbon.GetIdx(), nbr.GetIdx())
                    if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
                        is_carbonyl_neighbor = True
                        break
            if is_carbonyl_neighbor:
                break
        if is_carbonyl_neighbor:
            continue
        # If we get here, we found an extra hydroxyl.
        extra_oh_found = True
        break
    if not extra_oh_found:
        return False, "No extra hydroxy substituent found outside the carboxylic acid group"
    
    # If all tests pass, then the molecule can be considered a hydroxy fatty acid.
    reason = "Molecule contains a terminal carboxylic acid group, " 
    if n_carbons >= 8:
        reason += "a long aliphatic chain, "
    reason += "and extra hydroxy substituent(s), classifying it as a hydroxy fatty acid"
    return True, reason

# Example usage testing a few cases:
if __name__ == '__main__':
    test_smiles = [
        # True positives
        "C(C(O)=O)C[C@@H](/C=C/C=C\\C/C=C\\C/C=C\\C=C\\[C@H](C/C=C\\CC)O)O",  # resolvin D6
        "OC(=CC)C(O)=O",  # 2-hydroxy-2-butenoic acid
        "C(C(O)=O)C/C=C\\C[C@H](\\C=C\\C=C/C=C/C=C/[C@H]([C@@H](C/C=C\\CC)O)O)O",  # aspirin-triggered resolvin D2
        "CCCC\\C=C/C[C@@H]1O[C@H]1\\C=C\\[C@H](O)C\\C=C/CCCC(O)=O",  # (8R)-hydroxy-(11S,12S)-epoxyicosa-(5Z,9E,14Z)-trienoic acid
        "CC[C@@H](O)CCCCCCCCCCCCCCCC[C@@H](O)CC(O)=O",  # (3R,20R)-3,20-dihydroxyhenicosanoic acid
        "O[C@@H]([C@H](C)C(O)=O)C",  # 3-Hydroxy-2-methyl-[R-(R,S)]-butanoic acid
        "CCCCCC(O)CC(O)=O",  # 3-hydroxyoctanoic acid (should be accepted with relaxed aliphatic threshold)
        # False positive example (should be rejected as it has no extra –OH)
        "OC(=O)CCC(CCCC)CC",
    ]
    
    for s in test_smiles:
        res, reason = is_hydroxy_fatty_acid(s)
        print(f"SMILES: {s}\nResult: {res}, Reason: {reason}\n")