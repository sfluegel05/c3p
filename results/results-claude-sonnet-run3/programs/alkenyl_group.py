from rdkit import Chem
from rdkit.Chem import AllChem

def is_alkenyl_group(smiles: str):
    """
    Determines if a molecule is an alkenyl group (-CnH2n-1).
    
    Args:
        smiles (str): SMILES string of the molecule, with * representing attachment point
        
    Returns:
        bool: True if molecule is an alkenyl group, False otherwise
        str: Reason for classification
    """
    # Check for attachment point
    if '*' not in smiles:
        return False, "No attachment point (*) found"
        
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count atoms
    c_count = 0
    h_count = 0
    star_count = 0
    
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            c_count += 1
            h_count += len([n for n in atom.GetNeighbors() if n.GetSymbol() == 'H'])
            h_count += atom.GetNumImplicitHs()
        elif atom.GetSymbol() == '*':
            star_count += 1
            
    if star_count != 1:
        return False, "Must have exactly one attachment point"
        
    # Check for double bond
    double_bond_count = len([b for b in mol.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE])
    if double_bond_count == 0:
        return False, "No double bond found"
    elif double_bond_count > 1:
        return False, "Multiple double bonds found"
        
    # Check formula CnH(2n-1)
    expected_h = (2 * c_count - 1)
    if h_count != expected_h:
        return False, f"Incorrect H count. Expected {expected_h}, found {h_count}"
        
    # Check that attachment point is connected to carbon
    star_atom = None
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == '*':
            star_atom = atom
            break
            
    star_neighbors = [n for n in star_atom.GetNeighbors()]
    if len(star_neighbors) != 1 or star_neighbors[0].GetSymbol() != 'C':
        return False, "Attachment point must be connected to exactly one carbon"
        
    # All checks passed
    return True, f"Valid alkenyl group with formula C{c_count}H{h_count}"
# Pr=1.0
# Recall=0.8333333333333334