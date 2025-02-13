"""
Classifies: CHEBI:23899 icosanoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_icosanoid(smiles: str):
    """
    Determines if a molecule is an icosanoid based on its SMILES string.
    Icosanoids are signaling molecules derived from oxidation of C20 essential fatty acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an icosanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check the number of carbon atoms (20 carbon backbone)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (18 <= c_count <= 22):  # Allowing a small range to include derivatives
        return False, f"Carbon count ({c_count}) not consistent with C20 fatty acids"

    # Check for presence of oxygen atoms (indicative of oxidation products)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 2:
        return False, f"Insufficient number of oxygens ({o_count}) for an icosanoid"

    # Check for presence of double bonds (characteristic of EFAs)
    double_bonds = [bond for bond in mol.GetBonds() if bond.GetBondTypeAsDouble() == 2]
    if len(double_bonds) < 3:
        return False, f"Insufficient number of double bonds ({len(double_bonds)}) for an icosanoid"

    # Check for typical functional groups in icosanoids
    if not (Chem.MolFromSmarts("C=O") and Chem.MolFromSmarts("[OH]")).HasSubstructMatch(mol):
        return False, "Missing keto or hydroxyl groups typical of icosanoids"

    # Optional: Check for presence of cyclic structures like cyclopentane
    cyclic_rings = mol.GetRingInfo().NumRings()
    if cyclic_rings < 1:
        return False, "Missing cyclic structures typical in some icosanoids"
    
    return True, "Molecule fits the structural criteria of an icosanoid"