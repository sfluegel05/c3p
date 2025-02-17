"""
Classifies: CHEBI:24654 hydroxy fatty acid
"""
"""
Classifies: Hydroxy Fatty Acid
Definition: Any fatty acid carrying one or more hydroxy substituents.
A molecule is considered a hydroxy fatty acid if it:
  - Contains a terminal carboxylic acid group of the form C(=O)[OX2H] where the acid carbon is attached to exactly one other carbon.
  - Contains at least one extra –OH substituent that is not part of the acid group. Here, an extra –OH is defined as an O with one attached hydrogen that is bound to a carbon which is not carbonyl (not already involved in a C=O).
  - Has a minimum number of carbon atoms (we require at least 4 for any fatty acid).
  - For molecules with eight or more carbon atoms we additionally require that the fraction of heavy atoms (non‐H) that are carbon is at least 0.75 and that at least 85% of the carbon atoms are non‐aromatic.
  - Contains at most one ring (to allow for small epoxide rings but avoid polycyclic structures).
Molecules not meeting these criteria are not classified as hydroxy fatty acids.
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
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens (this will help in identifying -OH groups correctly)
    mol = Chem.AddHs(mol)
    
    # Allow at most one ring (small epoxides allowed)
    if mol.GetRingInfo().NumRings() > 1:
        return False, "Molecule has too many rings to be a fatty acid"
    
    # --- Find terminal carboxylic acid group ---
    # Use SMARTS for carboxylic acid: acid carbon (sp2) double bonded to an O and single bonded to an O with an H.
    acid_pattern = Chem.MolFromSmarts("[C;!R](=O)[O;H]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid group found (not a fatty acid)"
    
    valid_acid_match = None
    # In our pattern, match[0] is the acid carbon and match[1] is the hydroxyl oxygen.
    for match in acid_matches:
        acid_carbon = mol.GetAtomWithIdx(match[0])
        # Count carbon neighbours (should be exactly 1 for a terminal acid)
        carbon_neighbors = [nbr for nbr in acid_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) == 1:
            valid_acid_match = match
            break
    if valid_acid_match is None:
        return False, "No terminal carboxylic acid group found (acid carbon not attached to exactly one carbon)"
    acid_oh_idx = valid_acid_match[1]  # record acid hydroxyl index
    
    # --- Count carbon atoms in the molecule ---
    carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    n_carbons = len(carbons)
    if n_carbons < 4:
        return False, f"Too few carbons ({n_carbons}) to be considered a fatty acid"
    
    # --- For larger (more typical) fatty acids, require a high aliphatic fraction.
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() != 1]
    n_heavy = len(heavy_atoms)
    if n_carbons >= 8 and n_heavy > 0:  # apply only for longer chains
        n_heavy_carbon = sum(1 for atom in heavy_atoms if atom.GetAtomicNum() == 6)
        if (n_heavy_carbon / n_heavy) < 0.75:
            return False, "Molecule is not sufficiently aliphatic (low fraction of carbon among heavy atoms)"
        n_aliphatic_carbons = sum(1 for atom in carbons if not atom.GetIsAromatic())
        if (n_aliphatic_carbons / n_carbons) < 0.85:
            return False, "Molecule is not sufficiently aliphatic (too many aromatic carbons)"
    
    # --- Find extra hydroxy substituents (non-acid –OH) ---
    # Here we search for -OH groups that are not directly involved in a carbonyl.
    extra_oh_found = False
    for atom in mol.GetAtoms():
        # Look for oxygen atoms with exactly one bound hydrogen (as determined by explicit Hs)
        if atom.GetSymbol() == "O":
            # Check that it has an H (or count explicit hydrogens)
            H_count = sum(1 for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 1)
            if H_count != 1:
                continue
            # Exclude the acid -OH (by comparing indices)
            if atom.GetIdx() == acid_oh_idx:
                continue
            # Now require that the O is attached to a carbon that is not a carbonyl center.
            neighbors = atom.GetNeighbors()
            if not neighbors:
                continue
            attached = neighbors[0]
            if attached.GetAtomicNum() != 6:
                continue
            # Check if the attached carbon has a double-bonded oxygen (i.e. is part of a carbonyl).
            is_carbonyl = False
            for nbr in attached.GetNeighbors():
                bond = mol.GetBondBetweenAtoms(attached.GetIdx(), nbr.GetIdx())
                if nbr.GetAtomicNum() == 8 and bond is not None and bond.GetBondTypeAsDouble() == 1:
                    is_carbonyl = True
                    break
            if is_carbonyl:
                continue
            # Found an extra hydroxy not in the acid group.
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
        "OC(=CC)C(O)=O",  # 2-hydroxy-2-butenoic acid (4 carbons; accepted here despite its small size)
        "O[C@@H]([C@H](C)C(O)=O)C",  # 3-Hydroxy-2-methyl-[R-(R,S)]-butanoic acid (5 carbons)
        "C[C@@H](O)CCCCCCCCCCCCCCCC[C@@H](O)CC(O)=O",  # (3R,20R)-3,20-dihydroxyhenicosanoic acid
        # False positives:
        "OC(=O)CCC(CCCC)CC",  # 4-Ethyloctanoic acid: lacks any extra –OH; should be rejected.
    ]
    for s in test_smiles:
        res, reason = is_hydroxy_fatty_acid(s)
        print(f"SMILES: {s}\nResult: {res}, Reason: {reason}\n")