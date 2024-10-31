from rdkit import Chem
from rdkit.Chem import AllChem

def is_trimethoxyflavone(smiles: str):
    """
    Determines if a molecule is a trimethoxyflavone (flavone with exactly 3 methoxy groups).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a trimethoxyflavone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for flavone core structure (2-phenyl-4H-chromen-4-one)
    # More specific SMARTS pattern for flavone core
    flavone_pattern = Chem.MolFromSmarts('[#6]1=[#6]-c2c(o1)c(=O)cc1ccccc21')
    if not mol.HasSubstructMatch(flavone_pattern):
        return False, "No flavone core structure found"

    # Count methoxy groups (-OCH3)
    methoxy_pattern = Chem.MolFromSmarts('OC')
    methoxy_matches = mol.GetSubstructMatches(methoxy_pattern)
    
    # Count true methoxy groups
    true_methoxy_count = 0
    for match in methoxy_matches:
        o_atom = mol.GetAtomWithIdx(match[0])
        c_atom = mol.GetAtomWithIdx(match[1])
        
        # Check if carbon is CH3 (3 hydrogens and single bond only)
        if (c_atom.GetTotalNumHs() == 3 and 
            c_atom.GetDegree() == 1 and 
            o_atom.GetDegree() == 2):
            true_methoxy_count += 1

    if true_methoxy_count != 3:
        return False, f"Found {true_methoxy_count} methoxy groups, need exactly 3"

    return True, "Valid trimethoxyflavone with exactly 3 methoxy groups"
# Pr=None
# Recall=0.0