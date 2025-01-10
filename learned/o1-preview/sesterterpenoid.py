"""
Classifies: CHEBI:26660 sesterterpenoid
"""
"""
Classifies: CHEBI:26661 sesterterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is a sesterterpenoid based on its SMILES string.
    
    A sesterterpenoid is a terpenoid derived from a sesterterpene, which typically has 25 carbons
    built from five isoprene units (C5 units). The molecule may be rearranged or modified by 
    the removal of one or more skeletal atoms (generally methyl groups). Sesterterpenoids
    typically contain only carbon, hydrogen, and oxygen atoms.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a sesterterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for elements other than C, H, O
    allowed_atomic_nums = {1, 6, 8}  # H, C, O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Molecule contains heteroatom: {atom.GetSymbol()}"
    
    # Count number of carbons and hydrogens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    h_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
    
    # For sesterterpenoids, expect around 25 carbons, accept a small range due to modifications
    if c_count < 23 or c_count > 27:
        return False, f"Carbon count is {c_count}, which is not typical for a sesterterpenoid (expected 25±2 carbons)"
    
    # Estimate number of isoprene units
    isoprene_units = c_count / 5
    if not 4.5 <= isoprene_units <= 5.5:
        return False, f"Carbon count {c_count} does not correspond to 5 isoprene units"
    
    # Check for highly unsaturated structure (typical of terpenoids)
    unsaturation = rdMolDescriptors.CalcNumDoubleBonds(mol) + rdMolDescriptors.CalcNumAromaticRings(mol)
    if unsaturation < 2:
        return False, f"Low unsaturation ({unsaturation}); sesterterpenoids typically have multiple double bonds or rings"
    
    return True, "Molecule meets criteria for a sesterterpenoid"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:26661',
        'name': 'sesterterpenoid',
        'definition': 'Any terpenoid derived from a sesterterpene. The term includes compounds in which the C25 skeleton of the parent sesterterpene has been rearranged or modified by the removal of one or more skeletal atoms (generally methyl groups). Sometimes sesterterpenoids are erroneously referred to as sesterpenoids.',
        'parents': ['CHEBI:26873', 'CHEBI:26872']
    },
    'message': None,
    'success': True
}