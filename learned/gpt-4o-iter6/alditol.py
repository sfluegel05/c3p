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
    
    # Predominant carbon and oxygen atoms
    c_total = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    o_total = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'O')
    
    # Must have hydroxyl group for each central carbon
    if o_total < (c_total - 2): # excluding the two CH2OH termini
        return False, "Insufficient oxygen atoms for an alditol"

    # Enhanced check for terminal CH2OH groups (smarter pattern for ends)
    terminal_ch2oh_pattern = Chem.MolFromSmarts("[C;H2][O;X2H]")
    if len(mol.GetSubstructMatches(terminal_ch2oh_pattern)) < 2:
        return False, "Does not have typical terminal CH2OH groups"

    # Disallow other heteroatoms or complex structures that mimic polyols
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in ['C', 'O', 'H']:
            return False, f"Contains non CHO atoms, found: {atom.GetSymbol()}"

    return True, "Meets criteria: acyclic, terminal CH2OH groups and polyol characteristics found"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:17754',
        'name': 'alditol',
        'definition': "A carbohydrate that is an acyclic polyol having the general formula HOCH2[CH(OH)]nCH2OH, formally derivable from an aldose by reduction of the carbonyl group.",
        'parents': ['CHEBI:63221'],
    }
}