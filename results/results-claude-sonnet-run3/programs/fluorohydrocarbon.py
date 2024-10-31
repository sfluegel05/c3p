from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_fluorohydrocarbon(smiles: str):
    """
    Determines if a molecule is a fluorohydrocarbon (a hydrocarbon with one or more H atoms replaced by F).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fluorohydrocarbon, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Get all atoms
    atoms = mol.GetAtoms()
    
    # Track atom types present
    atom_symbols = set(atom.GetSymbol() for atom in atoms)
    
    # Must contain only C, H and F
    allowed_atoms = {'C', 'H', 'F'}
    if not atom_symbols.issubset(allowed_atoms):
        return False, f"Contains non C/H/F atoms: {atom_symbols - allowed_atoms}"
        
    # Must contain at least one C
    if 'C' not in atom_symbols:
        return False, "No carbon atoms present"
        
    # Must contain at least one F
    if 'F' not in atom_symbols:
        return False, "No fluorine atoms present"
        
    # Count number of F atoms
    f_count = sum(1 for atom in atoms if atom.GetSymbol() == 'F')
    
    # Count number of C atoms 
    c_count = sum(1 for atom in atoms if atom.GetSymbol() == 'C')
    
    # Verify carbon framework exists (carbons must be connected)
    c_atoms = [atom for atom in atoms if atom.GetSymbol() == 'C']
    visited = set()
    stack = [c_atoms[0].GetIdx()]
    
    while stack:
        current = stack.pop()
        if current not in visited:
            visited.add(current)
            atom = mol.GetAtomWithIdx(current)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'C':
                    stack.append(neighbor.GetIdx())
                    
    c_framework = len([idx for idx in visited if mol.GetAtomWithIdx(idx).GetSymbol() == 'C'])
    if c_framework != c_count:
        return False, "Carbon atoms are not connected in a single framework"

    return True, f"Fluorohydrocarbon with {f_count} fluorine atoms and {c_count} carbon atoms"
# Pr=1.0
# Recall=1.0