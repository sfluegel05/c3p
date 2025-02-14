"""
Classifies: CHEBI:60245 inositol phosphoceramide
"""
"""
Classifies: CHEBI:39179 inositol phosphoceramide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_inositol_phosphoceramide(smiles: str):
    """
    Determines if a molecule is an inositol phosphoceramide based on its SMILES string.
    An inositol phosphoceramide is a phosphosphingolipid with an inositol residue
    linked to a ceramide moiety via a phosphodiester bridge.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an inositol phosphoceramide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for inositol residue
    inositol_pattern = Chem.MolFromSmarts("[C@H]1[C@@H]([C@H]([C@@H]([C@@H]([C@H]1O)O)O)O)O")
    inositol_matches = mol.GetSubstructMatches(inositol_pattern)
    if len(inositol_matches) != 1:
        return False, "Exactly one inositol residue required"

    # Look for phosphodiester group
    phosphodiester_pattern = Chem.MolFromSmarts("O=P(OCC)(OC)")
    phosphodiester_matches = mol.GetSubstructMatches(phosphodiester_pattern)
    if len(phosphodiester_matches) != 1:
        return False, "Exactly one phosphodiester group required"

    # Look for ceramide moiety
    sphingoid_base_pattern = Chem.MolFromSmarts("[N;H2,H1]CCC[C@@H]([C@@H](C)CC)")
    fatty_acid_pattern = Chem.MolFromSmarts("CCCCCCCCCC(=O)")
    ceramide_pattern = sphingoid_base_pattern.GetMol().GetBonds() + fatty_acid_pattern.GetMol().GetBonds()
    ceramide_matches = mol.GetSubstructMatches(ceramide_pattern)
    if len(ceramide_matches) != 1:
        return False, "Exactly one ceramide moiety required"

    # Check molecular weight range for inositol phosphoceramides
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 600 or mol_wt > 1200:
        return False, "Molecular weight outside typical range for inositol phosphoceramides"

    return True, "Contains inositol residue linked to ceramide moiety via phosphodiester bridge"