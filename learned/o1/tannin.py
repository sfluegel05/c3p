"""
Classifies: CHEBI:26848 tannin
"""
"""
Classifies: CHEBI:33713 tannin
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tannin(smiles: str):
    """
    Determines if a molecule is a tannin based on its SMILES string.
    A tannin is defined as 'Any of a group of astringent polyphenolic vegetable principles or compounds,
    chiefly complex glucosides of catechol and pyrogallol.'
    This function checks for structural features characteristic of hydrolyzable and condensed tannins.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tannin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define patterns for hydrolyzable tannins (gallotannins and ellagitannins)
    # Galloyl group attached via ester linkage to sugar
    galloyl_ester_pattern = Chem.MolFromSmarts('O[C;R][C;R](O)COC(=O)c1cc(O)ccc1O')
    galloyl_ester_matches = mol.HasSubstructMatch(galloyl_ester_pattern)

    # Hexahydroxydiphenoyl (HHDP) group attached to sugar
    hhdp_ester_pattern = Chem.MolFromSmarts('O[C;R][C;R](O)COC(=O)c1c(O)cc2oc(=O)c(O)cc2c1O')
    hhdp_ester_matches = mol.HasSubstructMatch(hhdp_ester_pattern)

    # Define patterns for condensed tannins (proanthocyanidins)
    # Flavan-3-ol unit (catechin/epicatechin)
    flavan3ol_pattern = Chem.MolFromSmarts('[C@@H]1(O)Cc2cc(O)cc(O)c2O[C@@H]1c1cc(O)cc(O)c1')

    # Interflavan linkage C4-C8 or C4-C6
    c4_c8_linkage_pattern = Chem.MolFromSmarts('[C@@H]1(O)Cc2cc(O)cc(O)c2O[C@@H]1-c1c(O)cc(O)cc1')
    c4_c6_linkage_pattern = Chem.MolFromSmarts('[C@@H]1(O)Cc2cc(O)cc(O)c2O[C@@H]1-c1cc(O)cc(O)c1')

    flavan3ol_matches = mol.GetSubstructMatches(flavan3ol_pattern)
    c4_c8_matches = mol.GetSubstructMatches(c4_c8_linkage_pattern)
    c4_c6_matches = mol.GetSubstructMatches(c4_c6_linkage_pattern)

    # Check for hydrolyzable tannins
    if galloyl_ester_matches or hhdp_ester_matches:
        return True, "Molecule is a hydrolyzable tannin (gallotannin or ellagitannin)"

    # Check for condensed tannins
    if len(flavan3ol_matches) >= 2 and (len(c4_c8_matches) > 0 or len(c4_c6_matches) > 0):
        return True, "Molecule is a condensed tannin (proanthocyanidin) with multiple flavan-3-ol units linked via C-C bonds"

    # Additional check for larger condensed tannins
    if len(flavan3ol_matches) >= 3:
        return True, "Molecule is a larger condensed tannin with multiple flavan-3-ol units"

    # Check for multiple galloyl groups attached to a molecule
    galloyl_pattern = Chem.MolFromSmarts('c1cc(O)ccc1O[C](=O)')
    num_galloyl_groups = len(mol.GetSubstructMatches(galloyl_pattern))
    if num_galloyl_groups >= 2:
        return True, f"Molecule has multiple galloyl groups ({num_galloyl_groups}) attached via ester linkages"

    # If none of the above, molecule is likely not a tannin
    return False, "Molecule does not match structural features characteristic of tannins"