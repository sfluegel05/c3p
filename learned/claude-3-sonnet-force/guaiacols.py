"""
Classifies: CHEBI:134251 guaiacols
"""
"""
Classifies: CHEBI:51504 guaiacols
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_guaiacols(smiles: str):
    """
    Determines if a molecule is a guaiacol based on its SMILES string.
    A guaiacol is any phenol carrying an additional methoxy substituent at the ortho-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a guaiacol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for phenol ring pattern (aromatic ring with -OH group)
    phenol_pattern = Chem.MolFromSmarts("c1ccccc1O")
    if not mol.HasSubstructMatch(phenol_pattern):
        return False, "Not a phenol"
    
    # Find the phenolic -OH atom
    hydroxy_atom = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == "O" and atom.GetFormalCharge() == 0 and mol.GetBondBetweenAtoms(atom.GetIdx(), atom.GetNeighbors()[0]).IsInRingOfSize(6)][0]
    
    # Look for methoxy groups (-O-C) at ortho position relative to the phenolic -OH
    methoxy_pattern = Chem.MolFromSmarts("OC")
    methoxy_groups = mol.GetSubstructMatches(methoxy_pattern)
    
    has_ortho_methoxy = False
    for methoxy_match in methoxy_groups:
        methoxy_atom = methoxy_match[0]
        if mol.GetBondBetweenAtoms(methoxy_atom, hydroxy_atom).IsInRingOfSize(6):
            has_ortho_methoxy = True
            break
    
    if not has_ortho_methoxy:
        return False, "No ortho-methoxy group found relative to the phenolic -OH"
    
    return True, "Molecule contains a phenol with at least one ortho-methoxy substituent"