"""
Classifies: CHEBI:18133 hexose
"""
"""
Classifies: CHEBI:18133 hexose
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_hexose(smiles: str):
    """
    Determines if a molecule is a hexose based on its SMILES string.
    A hexose is any six-carbon monosaccharide which in its linear form contains either an aldehyde group at position 1 (aldohexose)
    or a ketone group at position 2 (ketohexose).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hexose, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Remove hydrogens for accurate atom counting
    mol = Chem.RemoveHs(mol)

    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 6:
        return False, f"Molecule has {c_count} carbon atoms, expected 6"

    # Count oxygen atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    # Hexoses typically have 6 carbons and 6 oxygens (C6H12O6)
    if o_count < 5 or o_count > 8:
        return False, f"Molecule has {o_count} oxygen atoms, expected between 5 and 8"

    # Check for glycosidic bonds (exclude molecules with glycosidic bonds)
    # Glycosidic bond pattern: anomeric carbon (connected to two oxygens) linked to another saccharide unit
    glycosidic_bond_pattern = Chem.MolFromSmarts("[C@H]1(O)[O][C@@H](O)[C@H](O)[C@@H](O)[C@H]1O")
    if mol.HasSubstructMatch(glycosidic_bond_pattern):
        return False, "Molecule appears to have glycosidic bonds, may not be a monosaccharide"

    # Check for aldehyde group (aldohexose)
    # Aldehyde pattern: [CHO]
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[#6]")
    has_aldehyde = mol.HasSubstructMatch(aldehyde_pattern)

    # Check for ketone group (ketohexose)
    # Ketone pattern: [#6][CX3](=O)[#6]
    ketone_pattern = Chem.MolFromSmarts("[#6][CX3](=O)[#6]")
    has_ketone = mol.HasSubstructMatch(ketone_pattern)

    # Check for cyclic hemiacetal (pyranose form)
    pyranose_pattern = Chem.MolFromSmarts("C1OC(O)C(O)C(O)C1O")
    has_pyranose = mol.HasSubstructMatch(pyranose_pattern)

    # Check for cyclic hemiketal (furanose form)
    furanose_pattern = Chem.MolFromSmarts("C1OC(O)C(O)C1O")
    has_furanose = mol.HasSubstructMatch(furanose_pattern)

    # Check for ring sizes (5 or 6-membered rings with oxygen)
    ring_info = mol.GetRingInfo()
    has_valid_ring = False
    for ring in ring_info.AtomRings():
        if len(ring) == 5 or len(ring) == 6:
            # Check if ring contains exactly one oxygen atom
            o_in_ring = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            if o_in_ring == 1:
                has_valid_ring = True
                break

    # Check for monosaccharide (no glycosidic bonds)
    # We can check that all oxygens are attached to carbons and not forming O-glycosidic linkages
    # Also, the molecule should not have phosphates or other non-standard substituents

    # Exclude molecules with phosphates, sulfates, or other large substituents
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in [1, 6, 8]:  # Allow only C, H, O
            return False, f"Molecule contains atom type {atom.GetSymbol()}, which is not typical in a hexose"

    # Evaluate if molecule has aldehyde, ketone, pyranose, or furanose form
    if has_aldehyde or has_ketone or has_pyranose or has_furanose or has_valid_ring:
        return True, "Molecule is a hexose: 6 carbons, appropriate oxygen content, and contains aldehyde/ketone or cyclic hemiacetal/hemiketal"
    else:
        return False, "Molecule does not contain aldehyde, ketone, or cyclic hemiacetal/hemiketal in expected positions"