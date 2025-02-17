"""
Classifies: CHEBI:24654 hydroxy fatty acid
"""
"""
Classifies: Hydroxy Fatty Acid
Definition: Any fatty acid carrying one or more hydroxy substituents.
A valid hydroxy fatty acid must have:
  - A terminal carboxylic acid group in the form C(=O)[OX2H] where the acid carbon is attached to exactly one other carbon.
  - At least one additional –OH group not part of the acid.
  - A sufficiently long aliphatic chain (we require at least 8 carbons).
  - A high fraction of heavy atoms (non‐H) are carbon and most of the carbon atoms (if any) are non‐aromatic.
  - Not too many rings (we allow at most one, which may occur in small epoxide groups).
Molecules not meeting these criteria will not be classified as hydroxy fatty acids.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a hydroxy fatty acid based on its SMILES string.

    Criteria:
      - Must contain a carboxylic acid group in the form C(=O)[OX2H].
      - The acid group must be terminal: the acid carbon must be attached to exactly one other carbon.
      - Must contain at least one additional –OH group (excluding the acid group OH).
      - The molecule must be mainly aliphatic:
           * It must have at least 8 carbon atoms.
           * The fraction of heavy (non-hydrogen) atoms that are carbon must be at least 0.75.
           * Additionally, the fraction of carbon atoms that are non‐aromatic must be at least 0.85.
      - The molecule should not have more than one ring.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a hydroxy fatty acid, False otherwise.
        str: Reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit H's to help in recognizing -OH groups.
    mol = Chem.AddHs(mol)
    
    # Check for rings: many fatty acids are acyclic (or nearly so).
    if mol.GetRingInfo().NumRings() > 1:
        return False, "Molecule has too many rings to be a fatty acid"
    
    # --- Check for terminal carboxylic acid group ---
    # SMARTS for carboxylic acid: the acid is represented as CO(=O)H 
    acid_pattern = Chem.MolFromSmarts("C(=O)[OX2H]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid group found (not a fatty acid)"
    
    # Filter acid matches to accept only terminal acid groups.
    valid_acid_match = None
    for match in acid_matches:
        # In our SMARTS "C(=O)[OX2H]", match[0] is the acid carbon.
        acid_carbon = mol.GetAtomWithIdx(match[0])
        # Count how many carbon neighbors are attached to the acid carbon.
        carbon_neighbors = [nbr for nbr in acid_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) == 1:
            valid_acid_match = match
            break
    if valid_acid_match is None:
        return False, "No terminal carboxylic acid group found (acid carbon not attached to exactly one carbon)"
    
    # Record the index of the acid hydroxyl (to exclude it from our extra –OH search).
    acid_oh_idx = valid_acid_match[1]  # second atom in the match is the acid hydroxyl oxygen

    # --- Check that the molecule has a sufficiently long aliphatic chain ---
    # Count all carbon atoms.
    carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    n_carbons = len(carbons)
    if n_carbons < 8:
        return False, f"Too few carbons ({n_carbons}) to be considered a fatty acid"

    # Check the fraction of heavy atoms (non-H) that are carbon.
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() != 1]
    n_heavy = len(heavy_atoms)
    n_heavy_carbon = sum(1 for atom in heavy_atoms if atom.GetAtomicNum() == 6)
    if n_heavy == 0 or (n_heavy_carbon / n_heavy) < 0.75:
        return False, "Molecule is not sufficiently aliphatic (low fraction of carbon among heavy atoms)"
    
    # Check that most carbon atoms are non-aromatic.
    n_aliphatic_carbons = sum(1 for atom in carbons if not atom.GetIsAromatic())
    if n_carbons > 0 and (n_aliphatic_carbons / n_carbons) < 0.85:
        return False, "Molecule is not sufficiently aliphatic (too many aromatic carbons)"
    
    # --- Check for an extra hydroxy substituent (outside the acid group) ---
    oh_pattern = Chem.MolFromSmarts("[OX2H]")
    oh_matches = mol.GetSubstructMatches(oh_pattern)
    
    extra_oh = []
    for match in oh_matches:
        # match is a tuple with the index of the O atom.
        if match[0] != acid_oh_idx:
            extra_oh.append(match[0])
    
    if not extra_oh:
        return False, "No hydroxy substituent found outside the carboxylic acid group"
    
    return True, "Molecule contains a terminal carboxylic acid group, a sufficiently long aliphatic chain, and extra hydroxy substituent(s), classifying it as a hydroxy fatty acid"


# Example usage:
if __name__ == '__main__':
    test_smiles = [
        # True positives (should be classified as hydroxy fatty acids)
        "C(C(O)=O)C[C@@H](/C=C/C=C\\C/C=C\\C/C=C\\C=C\\[C@H](C/C=C\\CC)O)O",  # resolvin D6
        "OC(=CC)C(O)=O",  # 2-hydroxy-2-butenoic acid  --> although small, we now require an 8-carbon chain so will be rejected
        "CCCCCCCC(O)C(O)=O",  # simple long chain acid with extra -OH (example not in list but helpful)
        # False positives / negatives
        "OC(=O)CCC(CCCC)CC",  # 4-Ethyloctanoic acid (lacks extra -OH outside the acid)
    ]
    for s in test_smiles:
        res, reason = is_hydroxy_fatty_acid(s)
        print(f"SMILES: {s}")
        print(f"Result: {res}, Reason: {reason}\n")