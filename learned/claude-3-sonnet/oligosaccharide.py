"""
Classifies: CHEBI:50699 oligosaccharide
"""
"""
Classifies: CHEBI:36973 oligosaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_oligosaccharide(smiles: str):
    """
    Determines if a molecule is an oligosaccharide based on its SMILES string.
    An oligosaccharide is a compound in which monosaccharide units are joined by glycosidic linkages,
    forming a defined structure as opposed to a polymer of unspecified length or a homologous mixture.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oligosaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for monosaccharides and glycosidic bonds
    monosaccharide_pattern = Chem.MolFromSmarts("[C@H]1[C@H]([C@@H]([C@H]([C@@H]1O)O)O)O")  # Hexose pattern
    glycosidic_bond_pattern = Chem.MolFromSmarts("[OX2][CR]")
    
    # Check for presence of multiple monosaccharide units
    monosaccharide_matches = mol.GetSubstructMatches(monosaccharide_pattern)
    if len(monosaccharide_matches) < 2:
        return False, "Fewer than 2 monosaccharide units found"
    
    # Check for glycosidic bonds connecting monosaccharide units
    glycosidic_bond_matches = mol.GetSubstructMatches(glycosidic_bond_pattern)
    if not glycosidic_bond_matches:
        return False, "No glycosidic bonds found"
    
    # Check if connected monosaccharide units form a defined structure
    # (This is a simplistic check based on molecular weight and number of components)
    n_components = Chem.Mol.GetNumAtomsSplitByConnectedComponent(mol)
    mol_wt = Chem.Descriptors.MolWt(mol)
    if n_components > 1 or mol_wt < 300:  # Arbitrary threshold for defined structure
        return False, "Likely a mixture or polymer, not a defined oligosaccharide structure"
    
    return True, "Contains multiple monosaccharide units joined by glycosidic linkages in a defined structure"