"""
Classifies: CHEBI:26935 tetraterpenoid
"""
"""
Classifies: CHEBI:26964 tetraterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tetraterpenoid(smiles: str):
    """
    Determines if a molecule is a tetraterpenoid based on its SMILES string.
    A tetraterpenoid is a terpenoid derived from a tetraterpene (C40 skeleton).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetraterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 35:
        return False, f"Too few carbons ({c_count}) for tetraterpenoid"
    if c_count > 45:
        return False, f"Too many carbons ({c_count}) for tetraterpenoid"

    # Count number of double bonds
    double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE)
    if double_bonds < 10:
        return False, f"Too few double bonds ({double_bonds}) for tetraterpenoid"

    # Check for presence of long conjugated system
    # Approximate by checking if the molecule has a high number of conjugated bonds
    conjugated_bonds = sum(1 for bond in mol.GetBonds() if bond.GetIsConjugated())
    if conjugated_bonds < 10:
        return False, f"Too few conjugated bonds ({conjugated_bonds}) for tetraterpenoid"

    # Check for terpenoid functional groups (presence of oxygen atoms)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count == 0:
        return False, "No oxygen atoms found; may not be a terpenoid"

    return True, "Molecule meets criteria for tetraterpenoid (C40 skeleton with terpenoid features)"