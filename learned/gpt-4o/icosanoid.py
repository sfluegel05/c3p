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

    # Check the number of carbon atoms (focus on broader range for C20 derivatives)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 18:  # Lower bound to ensure partial C20 structures or derivatives are considered
        return False, f"Carbon count ({c_count}) not sufficient for icosanoid derivatives"

    # Check for presence of oxygen atoms (indicative of oxidation patterns)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 3:
        return False, f"Insufficient number of oxygens ({o_count}), expect oxidized features in icosanoids"

    # Check for presence of double bonds (characteristic of EFAs)
    double_bonds = len([bond for bond in mol.GetBonds() if bond.GetBondTypeAsDouble() == 2])
    if double_bonds < 3:
        return False, f"Insufficient number of double bonds ({double_bonds}) for an icosanoid"

    # Check for typical functional groups: keto (C=O) and hydroxyl (OH)
    # Look for predominant icosanoid features rather than just presence
    keto_pattern = Chem.MolFromSmarts("C=O")
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    if not (mol.HasSubstructMatch(keto_pattern) and mol.HasSubstructMatch(hydroxyl_pattern)):
        return False, "Lacks both common functional groups (keto and hydroxyl) in majority of icosanoids"

    # Recognize cyclic elements (typical cyclopentane or structural analogs)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 0:
        return False, "Expected cyclic structure not detected in icosanoids"

    # Consider ester or amide linkages that are typical in complex icosanoids
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H0]")  # Ester linkage
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)")    # Amide linkage
    if mol.HasSubstructMatch(ester_pattern) or mol.HasSubstructMatch(amide_pattern):
        return True, "Structure contains substantially icosanoid-linked features such as ester or amide linkages"
    
    # Verify overall conformation align with EFAs conceptual ecosystem
    chemical_motif_alignment = Chem.MolFromSmarts("CCC=C")  # Simple motif patterns
    if not mol.HasSubstructMatch(chemical_motif_alignment):
        return False, "Molecule lacks aligned chemical motif patterns typical of icosanoids"

    return True, "Molecule generally fits refined structural criteria suggestive of an icosanoid"