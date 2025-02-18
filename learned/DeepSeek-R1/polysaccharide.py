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

    # Define monosaccharide pattern (pyranose/furanose with multiple hydroxyls)
    mono_pattern = Chem.MolFromSmarts("[#6]1(-[#8]-[#6]-[#6]-[#6]-[#6]-[#6]-1)-[OH]")
    matches = mol.GetSubstructMatches(mono_pattern)
    
    # Need at least 11 monosaccharide units (definition specifies >10)
    if len(matches) < 11:
        return False, f"Only {len(matches)} monosaccharide units found"

    # Check glycosidic linkages (O connecting two sugar rings)
    glycosidic_pattern = Chem.MolFromSmarts("[#6]-[#8]-[#6]1(-[#8]-[#6]-[#6]-[#6]-[#6]-[#6]-1)")
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No glycosidic linkages detected"

    # Check molecular weight (polysaccharides are typically >1500 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 1500:
        return False, f"Molecular weight ({mol_wt:.1f} Da) too low for polysaccharide"

    # Count oxygen atoms (approximate check for multiple glycosidic bonds)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o