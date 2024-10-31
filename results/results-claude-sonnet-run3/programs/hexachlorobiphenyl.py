from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_hexachlorobiphenyl(smiles: str):
    """
    Determines if a molecule is a hexachlorobiphenyl (C12H4Cl6).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a hexachlorobiphenyl, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Check molecular formula by counting atoms
    c_count = len([atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'C'])
    h_count = sum(atom.GetTotalNumHs() for atom in mol.GetAtoms())
    cl_count = len([atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'Cl'])
    
    if c_count != 12:
        return False, "Must contain exactly 12 carbon atoms"
    if h_count != 4:
        return False, "Must contain exactly 4 hydrogen atoms"
    if cl_count != 6:
        return False, "Must contain exactly 6 chlorine atoms"
    
    # Check for presence of biphenyl core structure
    biphenyl_pattern = Chem.MolFromSmarts('c1ccccc1-c1ccccc1')
    if not mol.HasSubstructMatch(biphenyl_pattern):
        return False, "Must contain biphenyl core structure"
        
    # Check that all carbons are aromatic
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'C']
    if not all(atom.GetIsAromatic() for atom in carbon_atoms):
        return False, "All carbons must be aromatic"
        
    # Check that chlorines are connected to aromatic carbons
    chlorine_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'Cl']
    for cl in chlorine_atoms:
        neighbors = cl.GetNeighbors()
        if len(neighbors) != 1 or not neighbors[0].GetIsAromatic() or neighbors[0].GetSymbol() != 'C':
            return False, "All chlorines must be connected to aromatic carbons"
            
    return True, "Valid hexachlorobiphenyl structure"
# Pr=1.0
# Recall=1.0