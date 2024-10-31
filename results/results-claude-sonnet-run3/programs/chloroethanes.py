from rdkit import Chem
from rdkit.Chem import AllChem

def is_chloroethanes(smiles: str):
    """
    Determines if a molecule is a chloroethane (ethane with one or more H replaced by Cl).
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a chloroethane, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Add hydrogens explicitly
    mol = Chem.AddHs(mol)
        
    # Count number of carbons
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    if num_carbons != 2:
        return False, "Not an ethane derivative - must have exactly 2 carbons"
        
    # Count number of chlorines
    num_chlorines = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'Cl')
    if num_chlorines == 0:
        return False, "No chlorine atoms present"
        
    # Check that molecule only contains C, H, Cl
    allowed_atoms = {'C', 'H', 'Cl'}
    actual_atoms = {atom.GetSymbol() for atom in mol.GetAtoms()}
    if not actual_atoms.issubset(allowed_atoms):
        return False, f"Contains disallowed atoms: {actual_atoms - allowed_atoms}"
        
    # Check that carbons are connected
    c_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == 'C']
    if len(c_atoms) == 2:
        path = Chem.GetShortestPath(mol, c_atoms[0], c_atoms[1])
        if len(path) != 2:  # Path length of 2 means directly connected
            return False, "Carbon atoms are not directly connected"
            
    # Verify total number of substituents
    total_bonds = sum(atom.GetTotalNumHs() + len(atom.GetBonds()) for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    if total_bonds != 8: # Each carbon should have 4 bonds
        return False, "Incorrect number of bonds on carbon atoms"
        
    return True, f"Chloroethane with {num_chlorines} chlorine atom(s)"
# Pr=1.0
# Recall=1.0