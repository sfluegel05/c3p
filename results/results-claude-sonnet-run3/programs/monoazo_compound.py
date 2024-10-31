from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdmolops import GetAdjacencyMatrix

def is_monoazo_compound(smiles: str):
    """
    Determines if a molecule is a monoazo compound (contains single -N=N- group)
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a monoazo compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Find all N=N bonds
    azo_bonds = []
    for bond in mol.GetBonds():
        atoms = [bond.GetBeginAtom(), bond.GetEndAtom()]
        if (bond.GetBondType() == Chem.BondType.DOUBLE and 
            all(atom.GetSymbol() == 'N' for atom in atoms)):
            azo_bonds.append(bond)
            
    # Check number of azo bonds
    if len(azo_bonds) == 0:
        return False, "No N=N bonds found"
    elif len(azo_bonds) > 1:
        return False, f"Multiple N=N bonds found ({len(azo_bonds)})"
        
    # For single azo bond, verify it's not part of a larger azoxy group
    azo_n_atoms = [azo_bonds[0].GetBeginAtom(), azo_bonds[0].GetEndAtom()]
    
    for n_atom in azo_n_atoms:
        for neighbor in n_atom.GetNeighbors():
            if neighbor.GetSymbol() == 'O' and neighbor.GetTotalNumHs() == 0:
                return False, "N=N bond is part of an azoxy group"
                
    return True, "Contains single N=N bond"
# Pr=None
# Recall=None