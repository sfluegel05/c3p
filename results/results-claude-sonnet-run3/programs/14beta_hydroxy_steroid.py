from rdkit import Chem
from rdkit.Chem import AllChem

def is_14beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 14beta-hydroxy steroid.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 14beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for steroid core (four fused rings)
    steroid_pattern = Chem.MolFromSmarts('[#6]1~[#6]~[#6]~[#6]2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]3~[#6]~2~[#6]~1')
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid core structure found"

    # Pattern for 14-beta-hydroxy steroid
    # The [C@@] specifies the beta stereochemistry at position 14
    # The pattern looks for the D ring with the beta-OH group at the junction
    beta_hydroxy_pattern = Chem.MolFromSmarts('[C]12[CH2][CH2][C]3[C]([C]1[CH2][CH2][C]1[C@@]([C]2)([C][C]3)O)[CH2][CH2][CH2][C]1')
    
    matches = mol.GetSubstructMatches(beta_hydroxy_pattern)
    if matches:
        # Verify the hydroxy group is present and has correct orientation
        for match in matches:
            # Get the carbon atom at position 14 (the one with OH group)
            c14_idx = match[8]  # Index of the carbon with OH in the SMARTS pattern
            c14_atom = mol.GetAtomWithIdx(c14_idx)
            
            # Check for OH group
            for neighbor in c14_atom.GetNeighbors():
                if neighbor.GetSymbol() == 'O' and neighbor.GetTotalNumHs() == 1:
                    return True, "Found 14-beta-hydroxy group in steroid structure"
    
    return False, "No 14-beta-hydroxy group found in correct position and orientation"
# Pr=None
# Recall=0.0