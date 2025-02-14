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
    This function checks for structural features characteristic of hydrolyzable and condensed tannins, 
    including galloyl groups attached to sugars, HHDP groups, and flavan-3-ol units linked via C-C or C-O bonds.

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

    # Define pattern for galloyl group attached via ester linkage to sugar (hydrolyzable tannins)
    galloyl_sugar_pattern = Chem.MolFromSmarts('O[C@H]1[*][C@H](O)[C@@H](O)[C@H](O)[C@@H]1O[C@H]1C(=O)Oc2cc(O)c(O)c(O)c2')

    # Define pattern for hexahydroxydiphenoyl (HHDP) group attached to sugar (ellagitannins)
    hhdp_sugar_pattern = Chem.MolFromSmarts('O[C@H]1[*][C@H](O)[C@@H](O)[C@@H](OC(=O)c2cc3oc(=O)c(O)cc3c(O)c2=O)[C@@H]1O')

    # Define pattern for flavan-3-ol units (catechin, epicatechin)
    flavan3ol_pattern = Chem.MolFromSmarts('Oc1cc(O)cc2OC[C@H]3OC(O)C[C@@H]3Oc12')

    # Define patterns for interflavan linkages in condensed tannins (proanthocyanidins)
    # C4-C8 linkage
    c4_c8_linkage_pattern = Chem.MolFromSmarts('Oc1ccc2O[C@@H]3COC[C@H]3Oc2c1-c1c(O)cc(O)cc1')
    # C4-C6 linkage
    c4_c6_linkage_pattern = Chem.MolFromSmarts('Oc1cc2O[C@@H]3COC[C@H]3Oc2c(O)c1-c1c(O)cc(O)cc1')

    # Define pattern for phloroglucinol units (phlorotannins)
    phloroglucinol_pattern = Chem.MolFromSmarts('Oc1cc(O)cc(O)c1')

    # Check for hydrolyzable tannins
    if mol.HasSubstructMatch(galloyl_sugar_pattern):
        return True, "Molecule is a gallotannin (hydrolyzable tannin with galloyl groups on sugar)"
    if mol.HasSubstructMatch(hhdp_sugar_pattern):
        return True, "Molecule is an ellagitannin (hydrolyzable tannin with HHDP group on sugar)"

    # Check for condensed tannins
    flavan3ol_matches = mol.GetSubstructMatches(flavan3ol_pattern)
    if len(flavan3ol_matches) >= 2:
        if mol.HasSubstructMatch(c4_c8_linkage_pattern) or mol.HasSubstructMatch(c4_c6_linkage_pattern):
            return True, "Molecule is a condensed tannin (proanthocyanidin) with flavan-3-ol units linked via interflavan bonds"

    # Check for phlorotannins (marine tannins composed of phloroglucinol units)
    phloroglucinol_matches = mol.GetSubstructMatches(phloroglucinol_pattern)
    if len(phloroglucinol_matches) >= 4:
        return True, "Molecule is a phlorotannin composed of multiple phloroglucinol units"

    # Check for catechin/epicatechin gallate esters
    catechin_gallate_pattern = Chem.MolFromSmarts('O[C@H]1[C@@H](O)C[C@@H](Oc2cc(O)cc(O)c12)OC(=O)c1cc(O)c(O)c(O)c1')
    if mol.HasSubstructMatch(catechin_gallate_pattern):
        return True, "Molecule is a catechin or epicatechin gallate ester"

    # Check for complex tannins (combination of hydrolyzable and condensed tannins)
    if (mol.HasSubstructMatch(galloyl_sugar_pattern) or mol.HasSubstructMatch(hhdp_sugar_pattern)) and len(flavan3ol_matches) >=1:
        return True, "Molecule is a complex tannin containing both hydrolyzable and condensed tannin units"

    # If none of the above, molecule is likely not a tannin
    return False, "Molecule does not match structural features characteristic of tannins"