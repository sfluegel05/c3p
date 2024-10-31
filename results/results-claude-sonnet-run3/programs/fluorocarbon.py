from rdkit import Chem
from rdkit.Chem import AllChem

def is_fluorocarbon(smiles: str):
    """
    Determines if a molecule is a fluorocarbon (contains only carbon and fluorine atoms).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a fluorocarbon, False otherwise
        str: Reason for classification
    """
    # Create RDKit molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Get all atoms in molecule
    atoms = mol.GetAtoms()
    
    # Check that only C and F atoms are present
    allowed_atoms = {'C', 'F'}
    atom_types = {atom.GetSymbol() for atom in atoms}
    
    if not atom_types.issubset(allowed_atoms):
        disallowed = atom_types - allowed_atoms
        return False, f"Contains disallowed atoms: {', '.join(disallowed)}"
        
    if 'F' not in atom_types:
        return False, "Contains no fluorine atoms"
        
    # Count number of each atom type
    c_count = sum(1 for atom in atoms if atom.GetSymbol() == 'C')
    f_count = sum(1 for atom in atoms if atom.GetSymbol() == 'F')
    
    return True, f"Fluorocarbon with {c_count} carbon atoms and {f_count} fluorine atoms"
# Pr=1.0
# Recall=1.0