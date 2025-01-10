"""
Classifies: CHEBI:26244 prenols
"""
"""
Classifies: prenols
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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
    
    # Check for hydroxyl group (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group (-OH) found"
    
    # Check for no rings (acyclic molecule)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings"
    
    # Check that molecule contains only carbon, hydrogen, and oxygen atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (1, 6, 8):
            return False, f"Contains atom other than C, H, or O: {atom.GetSymbol()}"
    
    # Count number of carbons
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Check that the number of carbons is at least 5 and a multiple of 5
    if num_carbons < 5:
        return False, f"Too few carbon atoms ({num_carbons}) to be a prenol"
    if num_carbons % 5 != 0:
        return False, f"Number of carbons ({num_carbons}) is not a multiple of 5"

    # Check for conjugated double bonds consistent with isoprene units
    num_double_bonds = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            num_double_bonds += 1
    min_expected_double_bonds = num_carbons // 5
    if num_double_bonds < min_expected_double_bonds:
        return False, f"Insufficient number of double bonds ({num_double_bonds}) for isoprene units"

    # Check if structure matches repeating isoprene units
    # Define SMARTS pattern for isoprene unit connected head-to-tail
    isoprene_pattern = Chem.MolFromSmarts("C(=C(C)C)C")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) < 1:
        return False, "No isoprene units found"
    
    return True, "Molecule is an acyclic alcohol with carbon skeleton composed of isoprene units"

__metadata__ = {   
    'chemical_class': {   
        'id': None,
        'name': 'prenols',
        'definition': "Any alcohol possessing the general formula H-[CH2C(Me)=CHCH2]nOH in which the carbon skeleton is composed of one or more isoprene units (biogenetic precursors of the isoprenoids).",
        'parents': []
    }
}