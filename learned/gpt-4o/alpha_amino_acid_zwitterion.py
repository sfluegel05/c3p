"""
Classifies: CHEBI:78608 alpha-amino-acid zwitterion
"""
from rdkit import Chem

def is_alpha_amino_acid_zwitterion(smiles: str):
    """
    Determines if a molecule is an alpha-amino-acid zwitterion based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is an alpha-amino-acid zwitterion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the presence of [NH3+] (protonated amine) and C([O-])=O (deprotonated carboxylate)
    amine_pattern = Chem.MolFromSmarts("[NH3+]")
    carboxylate_pattern = Chem.MolFromSmarts("C([O-])=O")

    if not mol.HasSubstructMatch(amine_pattern):
        return False, "No protonated amine group found"

    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No deprotonated carboxylate group found"
    
    # Check if amine and carboxylate are on the same carbon (alpha position)
    amine_matches = mol.GetSubstructMatches(amine_pattern)
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    
    alpha_match = False
    for amine_idx in amine_matches:
        for carboxylate_idx in carboxylate_matches:
            if any(mol.GetBondBetweenAtoms(amine_idx[0], carboxylate_idx[i]) for i in range(3)):
                alpha_match = True
                break
        if alpha_match:
            break

    if not alpha_match:
        return False, "Amine and carboxylate are not on the same carbon (no alpha configuration)"
    
    return True, "Contains alpha-amino acid structure with zwitterionic charges"

# Example usage:
# smiles = "CC(=O)CC([NH3+])C([O-])=O"  # Example SMILES from the list
# result, reason = is_alpha_amino_acid_zwitterion(smiles)
# print(f"Result: {result}, Reason: {reason}")