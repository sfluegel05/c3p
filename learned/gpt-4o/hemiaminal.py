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
    
    # Adjusted SMARTS pattern for hemiaminal: C with OH and NR groups
    hemiaminal_pattern = Chem.MolFromSmarts("[C;H2,H3;!R]([OH])[NH2,NHR,NHR2]")
    
    # Find matches with the hemiaminal pattern
    matches = mol.GetSubstructMatches(hemiaminal_pattern)
    if matches:
        for match in matches:
            carbon_index = match[0]
            oxygen_index = match[1]
            nitrogen_index = match[2]
            
            # Ensuring OH and NR are directly attached to the same carbon atom
            carbon_atom = mol.GetAtomWithIdx(carbon_index)
            if carbon_atom.GetDegree() != 3:
                continue  # Should have exactly 3 substituents (2 plus implicit hydrogen)
        
            # Ensure atoms are correctly classified
            oxygen_atom = mol.GetAtomWithIdx(oxygen_index)
            nitrogen_atom = mol.GetAtomWithIdx(nitrogen_index)
            if (oxygen_atom.GetAtomicNum() == 8 and 
                nitrogen_atom.GetAtomicNum() == 7 and
                carbon_atom.GetAtomicNum() == 6):
                return True, "Contains hemiaminal motif: C with OH and NR attached"
    
    return False, "No hemiaminal motif found"

# Example usage for debugging:
# result, reason = is_hemiaminal("NC(O)C(O)=O")
# print(f"Is hemiaminal: {result}, Reason: {reason}")