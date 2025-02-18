"""
Classifies: CHEBI:26848 tannin
"""
"""
Classifies: CHEBI:26848 tannin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tannin(smiles: str):
    """
    Determines if a molecule is a tannin based on its SMILES string.
    Tannins are polyphenolic compounds with multiple galloyl groups or phenolic subunits.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tannin, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for galloyl groups (O-C(=O)-attached trihydroxybenzene)
    galloyl_pattern = Chem.MolFromSmarts('[OX2]C(=O)c1c([OH])c([OH])c([OH])cc1')
    galloyl_count = len(mol.GetSubstructMatches(galloyl_pattern))
    if galloyl_count >= 2:
        return True, f"Contains {galloyl_count} galloyl groups"

    # Check for condensed tannin features (catechol/pyrogallol groups and high MW)
    catechol_pattern = Chem.MolFromSmarts('c1c([OH])c([OH])cccc1')  # Catechol (two adjacent OH)
    pyrogallol_pattern = Chem.MolFromSmarts('c1c([OH])c([OH])c([OH])ccc1')  # Pyrogallol (three adjacent OH)
    
    catechol_count = len(mol.GetSubstructMatches(catechol_pattern))
    pyrogallol_count = len(mol.GetSubstructMatches(pyrogallol_pattern))
    total_phenolic = catechol_count + pyrogallol_count

    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if total_phenolic >= 3 and mol_wt > 500:
        return True, f"Contains {total_phenolic} phenolic groups with MW {mol_wt:.1f}"

    # Check for flavan-3-ol like structures (common in condensed tannins)
    # This pattern approximates the flavan nucleus with multiple hydroxyls
    flavan_pattern = Chem.MolFromSmarts('[C@H]1OC2=C(C(=CC(=C2)O)O)C(C1)O')
    if mol.HasSubstructMatch(flavan_pattern):
        return True, "Contains flavan-3-ol subunit"

    return False, "Does not meet tannin criteria"