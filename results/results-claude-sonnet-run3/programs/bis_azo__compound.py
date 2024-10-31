from rdkit import Chem
from rdkit.Chem import AllChem

def is_bis_azo__compound(smiles: str):
    """
    Determines if a molecule is a bis(azo) compound containing two -N=N- groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bis(azo) compound, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Search for azo groups (-N=N-)
    # Using a SMARTS pattern to match azo groups
    azo_pattern = Chem.MolFromSmarts('[N]=[N]')
    matches = mol.GetSubstructMatches(azo_pattern)
    
    if len(matches) == 0:
        # Try alternative pattern for single bond representation
        azo_pattern_alt = Chem.MolFromSmarts('[N]-[N]')
        matches = mol.GetSubstructMatches(azo_pattern_alt)
        
    if len(matches) == 0:
        return False, "No azo groups found"
    elif len(matches) == 1:
        return False, "Only one azo group found"
    elif len(matches) >= 2:
        # Count valid azo groups (excluding those in azides or diazonium salts)
        valid_azo_count = 0
        for match in matches:
            n1, n2 = match
            atom1 = mol.GetAtomWithIdx(n1)
            atom2 = mol.GetAtomWithIdx(n2)
            
            # Check valence and formal charge
            if atom1.GetFormalCharge() == 0 and atom2.GetFormalCharge() == 0:
                if atom1.GetTotalValence() == 2 and atom2.GetTotalValence() == 2:
                    valid_azo_count += 1
        
        if valid_azo_count >= 2:
            return True, f"Contains {valid_azo_count} azo groups"
        else:
            return False, f"Contains {valid_azo_count} valid azo groups, need at least 2"
    
    return False, "No valid azo groups found"
# Pr=None
# Recall=None