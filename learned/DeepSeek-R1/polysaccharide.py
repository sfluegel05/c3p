"""
Classifies: CHEBI:18154 polysaccharide
"""
"""
Classifies: CHEBI:18154 polysaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polysaccharide(smiles: str):
    """
    Determines if a molecule is a polysaccharide based on its SMILES string.
    A polysaccharide is a biomacromolecule with >10 monosaccharide residues linked glycosidically.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polysaccharide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic sugar unit pattern: cyclic structure with multiple hydroxyls
    # Matches a 5 or 6-membered ring with at least 2 hydroxyl groups
    sugar_pattern = Chem.MolFromSmarts("[O][C@H]1[C@H](O)[C@H](O)[C@H](O)[C@H](O)O1")
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    
    # Check for at least 11 sugar units (definition specifies >10)
    if len(sugar_matches) < 11:
        return False, f"Only {len(sugar_matches)} sugar units found"

    # Glycosidic bond: oxygen connecting two anomeric carbons (C-O-C between rings)
    glycosidic_pattern = Chem.MolFromSmarts("[C@H]1[C@H](O)[C@H](O)[C@H](O)[C@H](O)O1-O-[C@H]2[C@H](O)[C@H](O)[C@H](O)[C@H](O)O2")
    glycosidic_bonds = mol.GetSubstructMatches(glycosidic_pattern)
    if len(glycosidic_bonds) < 10:  # At least 10 linkages for 11 units
        return False, f"Only {len(glycosidic_bonds)} glycosidic bonds detected"

    # Check molecular weight (approximate for >10 units: glucose polymer would be ~162*11=1782)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 1800:
        return False, f"Molecular weight ({mol_wt:.1f} Da) too low for polysaccharide"

    # Check oxygen count relative to number of units (each unit contributes multiple O)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    expected_o_min = 5 * len(sugar_matches)  # Approximate (each sugar has ~5 O)
    if o_count < expected_o_min * 0.8:  # Allow some variation
        return False, f"Oxygen count ({o_count}) inconsistent with polysaccharide structure"

    return True, f"Contains {len(sugar_matches)} sugar units with {len(glycosidic_bonds)} glycosidic bonds"