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
    
    # SMARTS pattern for hemiaminal: C bonded to OH and NRx (where R can vary)
    hemiaminal_pattern = Chem.MolFromSmarts("[C]([OH])([NH2,NHR,NHR2])")  
    
    # Check if the molecule matches the hemiaminal pattern
    matches = mol.GetSubstructMatches(hemiaminal_pattern)
    if matches:
        # Verify each match to ensure it is not part of larger functional groups
        for match in matches:
            carbon_atom = mol.GetAtomWithIdx(match[0])
            oxygen_atom = mol.GetAtomWithIdx(match[1])
            nitrogen_atom = mol.GetAtomWithIdx(match[2])
            
            # Ensure OH and NR group are both directly attached to the same carbon
            if carbon_atom.GetTotalValence() <= 4 and oxygen_atom.GetAtomicNum() == 8 and nitrogen_atom.GetAtomicNum() == 7:
                return True, "Contains hemiaminal motif: C with OH and NR attached"
    
    return False, "No hemiaminal motif found"

# Example usage:
# result, reason = is_hemiaminal("OC(N)CC")
# print(f"Is hemiaminal: {result}, Reason: {reason}")