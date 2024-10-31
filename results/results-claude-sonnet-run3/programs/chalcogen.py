from rdkit import Chem

def is_chalcogen(smiles: str):
    """
    Determines if a SMILES string represents a chalcogen atom (group 16 element).
    
    Args:
        smiles (str): SMILES string of the atom/molecule
        
    Returns:
        bool: True if atom is a chalcogen, False otherwise
        str: Reason for classification
    """
    # List of chalcogen elements (group 16)
    chalcogens = {'O', 'S', 'Se', 'Te', 'Po', 'Lv'}
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Check if molecule contains exactly one atom
    if mol.GetNumAtoms() != 1:
        return False, "Input must be a single atom"
        
    # Get the atom symbol
    atom = mol.GetAtomWithIdx(0)
    symbol = atom.GetSymbol()
    
    # Check if atom is a chalcogen
    if symbol in chalcogens:
        return True, f"{symbol} is a chalcogen (group 16 element)"
    else:
        return False, f"{symbol} is not a chalcogen"
# Pr=1.0
# Recall=1.0