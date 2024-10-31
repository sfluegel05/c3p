from rdkit import Chem
from rdkit.Chem import AllChem

def is_acyl_chloride(smiles: str):
    """
    Determines if a molecule is an acyl chloride (compound with R-C(=O)-Cl group).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an acyl chloride, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # SMARTS pattern for acyl chloride group: R-C(=O)-Cl
    pattern = Chem.MolFromSmarts('[#6](=O)Cl')
    
    if mol.HasSubstructMatch(pattern):
        matches = mol.GetSubstructMatches(pattern)
        num_matches = len(matches)
        
        # Get the carbon atom that's part of C(=O)Cl group
        c_atoms = [match[0] for match in matches]
        
        # Get the number of acyl chloride groups
        if num_matches == 1:
            return True, "Contains one acyl chloride group"
        else:
            return True, f"Contains {num_matches} acyl chloride groups"
            
    return False, "No acyl chloride group found"
# Pr=1.0
# Recall=1.0