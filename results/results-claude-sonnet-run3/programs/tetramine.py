from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tetramine(smiles: str):
    """
    Determines if a molecule is a tetramine (contains exactly 4 amino groups)
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a tetramine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
        
    # Find all nitrogen atoms
    nitrogen_atoms = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N':
            nitrogen_atoms.append(atom)
            
    if len(nitrogen_atoms) < 4:
        return False, "Contains fewer than 4 nitrogen atoms"
        
    # Count amino groups (-NH2, -NH-, -N<)
    amino_count = 0
    for n_atom in nitrogen_atoms:
        # Get number of hydrogens (explicit + implicit)
        h_count = n_atom.GetTotalNumHs()
        
        # Get number of non-H neighbors
        non_h_neighbors = len([a for a in n_atom.GetNeighbors()])
        
        # Amino group criteria:
        # -NH2: 2H, 1 non-H neighbor
        # -NH-: 1H, 2 non-H neighbors  
        # -N<: 0H, 3 non-H neighbors
        if ((h_count == 2 and non_h_neighbors == 1) or
            (h_count == 1 and non_h_neighbors == 2) or
            (h_count == 0 and non_h_neighbors == 3)):
            amino_count += 1
            
    if amino_count == 4:
        return True, "Contains exactly 4 amino groups"
    elif amino_count < 4:
        return False, f"Contains only {amino_count} amino groups"
    else:
        return False, f"Contains {amino_count} amino groups (more than 4)"
# Pr=0.8125
# Recall=1.0