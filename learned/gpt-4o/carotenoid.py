"""
Classifies: CHEBI:23044 carotenoid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

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
    if c_count < 30 or c_count > 50:  # Allow some flexibility around 40 carbons
        return False, f"Expected ~40 carbons, but found {c_count}"

    # Check for conjugated system (presence of alternating single and double bonds)
    max_conj_bonds = 0
    bonds = mol.GetBonds()
    current_conj_bonds = 0
    for bond in bonds:
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            current_conj_bonds += 1
        elif bond.GetBondType() == Chem.rdchem.BondType.SINGLE and current_conj_bonds > 0:
            current_conj_bonds += 1
        else:
            max_conj_bonds = max(max_conj_bonds, current_conj_bonds)
            current_conj_bonds = 0
    max_conj_bonds = max(max_conj_bonds, current_conj_bonds)

    if max_conj_bonds < 10:  # Require a minimum number of conjugated bonds
        return False, "Insufficient conjugated double bonds"

    # Check for potential modifications (OH, =O, epoxide)
    num_oxygen = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if num_oxygen > 10:  # Possible oxidation or hydroxylation, but modifiable
        return False, "Too many oxygens - might not be a carotenoid"
    
    # Look for cyclic modifications or long conjugated linear chains
    rings = mol.GetRingInfo().AtomRings()
    if len(rings) < 2:  # Often will involve small rings in additions
        return False, "Insufficient cyclic modifications"

    return True, "Matches the structural criteria of carotenoids"