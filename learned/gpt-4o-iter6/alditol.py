"""
Classifies: CHEBI:17522 alditol
"""
from rdkit import Chem

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

    # Checking for linear hydroxylated backbone
    # Ensuring primarily acyclic with small allowable rings
    c_total = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    o_total = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'O')
    
    # Minimal conditions for being a polyol (lots of OH groups in relation to C)
    if o_total < c_total / 2:
        return False, "Insufficient oxygen atoms for an alditol"

    # Must not contain disaccharide-like patterns heavily repeated
    glyco_pattern = Chem.MolFromSmarts("OC[C@H]1O[C@@H](CO)[C@@H](O)[C@H](O)[C@H]1O")  
    if mol.HasSubstructMatch(glyco_pattern):
        return False, "Contains repeated glycosidic-like units, not typical of simple alditols"

    # Ensure terminal CH2OH ends...
    terminal_oh_pattern = Chem.MolFromSmarts("CO")
    terminal_oh_matches = mol.GetSubstructMatches(terminal_oh_pattern)
    if len(terminal_oh_matches) < 2:
        return False, "Does not have typical terminal CH2OH groups"

    # Avoid carbonyl presence (identity thief of a hydride - alditol)
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=O")
    if mol.HasSubstructMatch(carbonyl_pattern):
        return False, "Contains carbonyl groups, indication of precursor not alditol"
    
    return True, "Meets criteria: terminal CH2OH groups and polyol characteristics found"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:17754',
        'name': 'alditol',
        'definition': "A carbohydrate that is an acyclic polyol having the general formula HOCH2[CH(OH)]nCH2OH, formally derivable from an aldose by reduction of the carbonyl group.",
        'parents': ['CHEBI:63221'],
    }
}