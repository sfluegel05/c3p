"""
Classifies: CHEBI:61778 triterpenoid saponin
"""
"""
Classifies: CHEBI:75214 triterpenoid saponin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_triterpenoid_saponin(smiles: str):
    """
    Determines if a molecule is a triterpenoid saponin based on its SMILES string.
    A triterpenoid saponin is a terpene glycoside with a triterpenoid backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triterpenoid saponin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for triterpenoid backbone (30 carbons, typically with multiple rings)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 30:
        return False, "Not enough carbons for a triterpenoid backbone"

    # Check for glycosidic linkages (presence of sugar moieties)
    glycosidic_pattern = Chem.MolFromSmarts("[C;!$(C=O)]-[O;!$(O=*)]-[C;!$(C=O)]")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)
    if len(glycosidic_matches) == 0:
        return False, "No glycosidic linkages found"

    # Check for multiple rings (typical of triterpenoids)
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    if n_rings < 4:
        return False, "Not enough rings for a triterpenoid backbone"

    # Check for hydroxyl groups (common in triterpenoids and saponins)
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 2:
        return False, "Not enough hydroxyl groups for a saponin"

    # Check for carboxylic acid groups (common in triterpenoids)
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if len(carboxyl_matches) == 0:
        return False, "No carboxylic acid groups found"

    # Check molecular weight (triterpenoid saponins are typically large)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for a triterpenoid saponin"

    return True, "Contains a triterpenoid backbone with glycosidic linkages"