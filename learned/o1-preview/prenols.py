"""
Classifies: CHEBI:26244 prenols
"""
"""
Classifies: prenols
"""
from rdkit import Chem

def is_prenols(smiles: str):
    """
    Determines if a molecule is a prenol based on its SMILES string.
    A prenol is any alcohol possessing the general formula H-[CH2C(Me)=CHCH2]nOH
    in which the carbon skeleton is composed of one or more isoprene units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prenol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for acyclic structure
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings"
    
    # Check that molecule contains only carbon, hydrogen, and oxygen atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (1, 6, 8):
            return False, f"Contains atom other than C, H, or O: {atom.GetSymbol()}"
    
    # Check that all bonds to oxygen are single bonds (exclude carbonyls)
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if begin_atom.GetAtomicNum() == 8 or end_atom.GetAtomicNum() == 8:
                return False, "Contains carbonyl group (C=O)"
    
    # Check for hydroxyl group (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group (-OH) found"
    
    # Check for primary alcohol group (-CH2OH)
    primary_alcohol_pattern = Chem.MolFromSmarts("[CX4H2][OX2H]")
    if not mol.HasSubstructMatch(primary_alcohol_pattern):
        return False, "No primary alcohol group (-CH2OH) found"
    
    # Check that molecule is acyclic hydrocarbon chain with hydroxyl group(s)
    # and does not contain other common functional groups (aldehydes, ketones, carboxylic acids)
    # We have already excluded molecules with C=O groups above
    
    # All conditions satisfied
    return True, "Molecule is an acyclic alcohol consistent with prenol structure"

__metadata__ = {   
    'chemical_class': {   
        'id': None,
        'name': 'prenols',
        'definition': "Any alcohol possessing the general formula H-[CH2C(Me)=CHCH2]nOH in which the carbon skeleton is composed of one or more isoprene units (biogenetic precursors of the isoprenoids).",
        'parents': []
    }
}