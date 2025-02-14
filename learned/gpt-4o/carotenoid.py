"""
Classifies: CHEBI:23044 carotenoid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors, rdchem

def is_carotenoid(smiles: str):
    """
    Determines if a molecule is a carotenoid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a carotenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 30 or c_count > 50:
        return False, f"Expected ~40 carbons, found {c_count}"

    # Look for conjugation - counted through alternating single and double bonds
    max_conj_bonds = 0
    bonds = mol.GetBonds()
    current_conj_bonds = 0
    for bond in bonds:
        if bond.GetBondType() == rdchem.BondType.DOUBLE:
            current_conj_bonds += 1
        elif bond.GetBondType() == rdchem.BondType.SINGLE and current_conj_bonds > 0:
            current_conj_bonds += 1
        else:
            max_conj_bonds = max(max_conj_bonds, current_conj_bonds)
            current_conj_bonds = 0
    max_conj_bonds = max(max_conj_bonds, current_conj_bonds)

    if max_conj_bonds < 8:
        return False, "Insufficient conjugated double bonds"

    # Check for oxygen atom presence (functional diversity with oligoethers or keto groups)
    num_oxygen = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if num_oxygen > 12:
        return False, "Unlikely carotenoid due to excess oxygen"

    # Assess for ring presence or equivalent structural complexity
    ring_info = mol.GetRingInfo().AtomRings()
    if len(ring_info) < 1:
        return False, "No cyclic structure found, generally rare"

    # Check for presence of potential carotenoid functional groups (OH, =O, epoxide)
    hydroxyl_group = Chem.MolFromSmarts('[OH]')
    if not mol.HasSubstructMatch(hydroxyl_group):
        return False, "No hydroxyl groups detected, uncommon in xanthophylls"

    return True, "Matches the structural criteria of carotenoids"