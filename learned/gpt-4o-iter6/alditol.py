"""
Classifies: CHEBI:17522 alditol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alditol(smiles: str):
    """
    Determines if a molecule is an alditol based on its SMILES string.
    An alditol is an acyclic polyol with the general formula 
    HOCH2[CH(OH)]nCH2OH (formally derived from an aldose by reduction).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alditol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure no cyclic structures
    if rdMolDescriptors.CalcNumRings(mol) > 0:
        return False, "Contains cyclic structures, alditol must be acyclic"
    
    # Count the carbon and oxygen atoms
    c_total = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    o_total = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'O')
    
    # Minimal condition for being a polyol: multiple oxygen (OH) per carbon
    if o_total < c_total / 2:
        return False, "Insufficient oxygen atoms for an alditol"

    # Enhanced check for terminal CH2OH groups
    terminal_ch2oh_pattern = Chem.MolFromSmarts("[CH2][OH]")
    if len(mol.GetSubstructMatches(terminal_ch2oh_pattern)) < 2:
        return False, "Does not have typical terminal CH2OH groups"

    # Exclude carbonyl presence to confirm alditol nature
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=O")
    if mol.HasSubstructMatch(carbonyl_pattern):
        return False, "Contains carbonyl groups, indication of precursor not alditol"

    # Additional check: avoid extensive glycosidic linkages with complex polyol
    glyco_pattern = Chem.MolFromSmarts("OC[C@H]1O[C@@H](CO)[C@@H](O)[C@@H](O)[C@H]1O")
    if mol.HasSubstructMatch(glyco_pattern):
        return False, "Contains repeated glycosidic-like units, not typical of simple alditols"
    
    return True, "Meets criteria: acyclic, terminal CH2OH groups and polyol characteristics found"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:17754',
        'name': 'alditol',
        'definition': "A carbohydrate that is an acyclic polyol having the general formula HOCH2[CH(OH)]nCH2OH, formally derivable from an aldose by reduction of the carbonyl group.",
        'parents': ['CHEBI:63221'],
    }
}