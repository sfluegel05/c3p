"""
Classifies: CHEBI:73080 hemiaminal
"""
from rdkit import Chem

def is_hemiaminal(smiles: str):
    """
    Determines if a molecule is a hemiaminal based on its SMILES string.
    A hemiaminal has an amino group and a hydroxy group attached to the same carbon atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hemiaminal, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"
    
    # Refined SMARTS pattern for hemiaminal (C bonded to any N and O atom)
    hemiaminal_pattern = Chem.MolFromSmarts("[#6][OH][#7]")  # More flexible, catering for any nitrogen attachment
    
    # Check if the molecule matches the hemiaminal pattern
    matches = mol.GetSubstructMatches(hemiaminal_pattern)
    if matches:
        # Further validate matches - ensure C-N and C-O are not part of different functionalities
        for match in matches:
            carbon_atom = mol.GetAtomWithIdx(match[0])
            oxygen_atom = mol.GetAtomWithIdx(match[1])
            nitrogen_atom = mol.GetAtomWithIdx(match[2])
            
            # Ensure OH and NHx are both attached to the same central C atom
            if carbon_atom.GetTotalNumHs() >= 1 and nitrogen_atom.GetTotalNumHs() >= 1:
                return True, "Contains a hemiaminal motif: C with OH and NHx attached"
    
    return False, "No hemiaminal motif found"

# Example usage:
# result, reason = is_hemiaminal("OC(N)CC")
# print(f"Is hemiaminal: {result}, Reason: {reason}")