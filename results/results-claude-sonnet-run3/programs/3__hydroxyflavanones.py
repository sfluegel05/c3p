from rdkit import Chem
from rdkit.Chem import AllChem

def is_3__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 3'-hydroxyflavanone (has a hydroxy substituent at position 3' of the phenyl ring)
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 3'-hydroxyflavanone, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for basic flavanone core structure
    flavanone_pattern = Chem.MolFromSmarts('[O;R1]1[CH2,CH1]C(=O)c2[cR1][cR1][cR1][cR1][cR1]2[CH1]1[c;R0]1[cR1][cR1][cR1][cR1][cR1][cR1]1')
    
    if not mol.HasSubstructMatch(flavanone_pattern):
        return False, "Not a flavanone structure"

    # Get the matches for the flavanone core
    matches = mol.GetSubstructMatches(flavanone_pattern)
    
    for match in matches:
        # The phenyl ring carbons will be the last 6 atoms in the match
        phenyl_ring = match[-6:]
        
        # Check each carbon in phenyl ring for hydroxy substituent
        for carbon_idx in phenyl_ring:
            carbon = mol.GetAtomWithIdx(carbon_idx)
            
            # Check neighbors for oxygen
            for neighbor in carbon.GetNeighbors():
                if neighbor.GetSymbol() == 'O' and neighbor.GetTotalNumHs() == 1:
                    # Check if this is the meta position (3')
                    if carbon_idx == phenyl_ring[2]: # Index 2 corresponds to 3' position
                        return True, "Found hydroxy group at 3' position"

    return False, "No hydroxy group found at 3' position"
# Pr=None
# Recall=0.0