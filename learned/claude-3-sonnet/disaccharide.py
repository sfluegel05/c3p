"""
Classifies: CHEBI:36233 disaccharide
"""
"""
Classifies: CHEBI:36973 disaccharide

A disaccharide is a compound in which two monosaccharides are joined by a glycosidic bond.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

# SMARTS patterns for common monosaccharides
MONOSACCHARIDE_PATTERNS = [
    "[OX2][CX4](O)[CX3](O)[CX3](O)[CX3](O)[CX3](O)[CX2]",  # Linear form
    "[OX2][CX4]1(O)[CX3]([CH2])[CX3]([CH2])[CX3](O)[CX3](O)[CX3]1",  # Pyranose form
    "[OX2][CX4]1(O)[CX3]([CH2])[CX3](O)[CX3](O)[CX3]1"  # Furanose form
]

def is_disaccharide(smiles: str) -> tuple[bool, str]:
    """
    Determines if a molecule is a disaccharide based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a disaccharide, False otherwise
        str: Reason for the classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule contains two monosaccharide subunits
    monosaccharide_units = []
    for pattern in MONOSACCHARIDE_PATTERNS:
        matches = mol.GetSubstructMatches(Chem.MolFromSmarts(pattern))
        monosaccharide_units.extend(matches)

    if len(monosaccharide_units) != 2:
        return False, f"Found {len(monosaccharide_units)} monosaccharide units, expected 2"

    # Check for a glycosidic bond connecting the two monosaccharide units
    glycosidic_bond_pattern = Chem.MolFromSmarts("[OX2][CX4]")
    bond_matches = mol.GetSubstructMatches(glycosidic_bond_pattern)

    # Check if the glycosidic bond connects the two monosaccharide units
    for bond_match in bond_matches:
        bond = mol.GetBondBetweenAtoms(bond_match[0], bond_match[1])
        if bond.GetBeginAtomIdx() in monosaccharide_units and bond.GetEndAtomIdx() in monosaccharide_units:
            break
    else:
        return False, "No glycosidic bond found between monosaccharide units"

    # Check molecular weight range (300-600 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 600:
        return False, f"Molecular weight ({mol_wt:.2f} Da) out of the expected range for disaccharides"

    return True, "Contains two monosaccharide units connected by a glycosidic bond"