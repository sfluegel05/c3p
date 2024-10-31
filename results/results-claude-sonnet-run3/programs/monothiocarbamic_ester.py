from rdkit import Chem
from rdkit.Chem import AllChem

def is_monothiocarbamic_ester(smiles: str):
    """
    Determines if a molecule is a monothiocarbamic ester - a thiocarbamic ester formally derived from a monothiocarbamic acid.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a monothiocarbamic ester, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Look for C(=O)S or C(=S)O pattern
    pattern1 = Chem.MolFromSmarts('[CX3](=[OX1])[SX2]') # C(=O)S pattern
    pattern2 = Chem.MolFromSmarts('[CX3](=[SX1])[OX2]') # C(=S)O pattern
    
    matches1 = mol.GetSubstructMatches(pattern1)
    matches2 = mol.GetSubstructMatches(pattern2)
    
    if not (matches1 or matches2):
        return False, "No C(=O)S or C(=S)O group found"

    # Check for nitrogen attached to the carbonyl/thiocarbonyl carbon
    if matches1:
        c_idx = matches1[0][0]
        pattern_n = Chem.MolFromSmarts('[NX3][CX3](=[OX1])[SX2]')
    else:
        c_idx = matches2[0][0]
        pattern_n = Chem.MolFromSmarts('[NX3][CX3](=[SX1])[OX2]')
        
    matches_n = mol.GetSubstructMatches(pattern_n)
    
    if not matches_n:
        return False, "No nitrogen atom attached to C(=O)S or C(=S)O group"

    # Get the central carbon atom and check its neighbors
    c_atom = mol.GetAtomWithIdx(c_idx)
    
    # Count number of double bonds
    double_bond_count = sum(1 for bond in c_atom.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE)
    
    if double_bond_count != 1:
        return False, "Central carbon should have exactly one double bond"

    if matches1:
        return True, "Contains C(=O)S-R group with attached nitrogen"
    else:
        return True, "Contains C(=S)O-R group with attached nitrogen"
# Pr=1.0
# Recall=1.0