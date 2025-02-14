"""
Classifies: CHEBI:26660 sesterterpenoid
"""
"""
Classifies: CHEBI:33551 sesterterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is a sesterterpenoid based on its SMILES string.
    A sesterterpenoid is a terpenoid with a C25 backbone, derived from a sesterterpene.

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

    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 22 or c_count > 28:
        return False, f"Molecule has {c_count} carbon atoms, outside the expected range for sesterterpenoids (22-28)"

    # Check for terpenoid-like structure
    try:
        Chem.SanitizeMol(mol)
        Chem.AssignStereochemistry(mol)
    except ValueError:
        return False, "Could not assign stereochemistry"

    # Check for rings and double bonds
    n_rings = mol.GetRingInfo().NumRings()
    n_double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)

    # Sesterterpenoids typically have 4-7 rings and 4-8 double bonds
    if n_rings < 3 or n_rings > 8:
        return False, f"Unexpected number of rings ({n_rings}), sesterterpenoids typically have 3-8 rings"
    if n_double_bonds < 3 or n_double_bonds > 10:
        return False, f"Unexpected number of double bonds ({n_double_bonds}), sesterterpenoids typically have 3-10 double bonds"

    # Check molecular weight - sesterterpenoids typically 300-600 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 600:
        return False, f"Molecular weight ({mol_wt:.2f} Da) out of typical sesterterpenoid range (300-600 Da)"

    # Additional check for rearranged or modified sesterterpene backbone
    sesterterpene_pattern = Chem.MolFromSmarts("[C&Ring1;R1][C&Ring1;R1][C&Ring1;R1][C&Ring1;R1][C&Ring1;R1][C&Ring1;R1][C&Ring1;R1][C&Ring1;R1][C&Ring1;R1][C&Ring1;R1][C&Ring1;R1][C&Ring1;R1][C&Ring1;R1][C&Ring1;R1][C&Ring1;R1][C&Ring1;R1][C&Ring1;R1][C&Ring1;R1][C&Ring1;R1][C&Ring1;R1][C&Ring1;R1][C&Ring1;R1][C&Ring1;R1][C&Ring1;R1][C&Ring1;R1]")
    if not mol.HasSubstructMatch(sesterterpene_pattern):
        return False, "No sesterterpene-like backbone found"

    return True, "Molecule matches structural characteristics of a sesterterpenoid"