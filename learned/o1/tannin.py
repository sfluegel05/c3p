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
    This function checks for the presence of multiple catechol or pyrogallol units, phenolic hydroxyl groups,
    and possible presence of sugar moieties connected via glycosidic or ester bonds.

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

    # Check for phenolic hydroxyl groups attached to aromatic rings
    phenol_pattern = Chem.MolFromSmarts('a O')
    phenol_matches = mol.GetSubstructMatches(phenol_pattern)
    num_phenol_groups = len(phenol_matches)

    # Tannins typically have multiple phenolic hydroxyl groups
    if num_phenol_groups < 6:
        return False, f"Only {num_phenol_groups} phenolic hydroxyl groups found, need at least 6"

    # Check for catechol units (benzene ring with two adjacent hydroxyl groups)
    catechol_pattern = Chem.MolFromSmarts('c1c(O)cc(O)cc1')
    num_catechol_units = len(mol.GetSubstructMatches(catechol_pattern))

    # Check for pyrogallol units (benzene ring with three adjacent hydroxyl groups)
    pyrogallol_pattern = Chem.MolFromSmarts('c1c(O)cc(O)cc1O')
    num_pyrogallol_units = len(mol.GetSubstructMatches(pyrogallol_pattern))

    total_catechol_pyrogallol = num_catechol_units + num_pyrogallol_units

    if total_catechol_pyrogallol < 2:
        return False, f"Only {total_catechol_pyrogallol} catechol/pyrogallol units found, need at least 2"

    # Check for sugar moieties (e.g., glucose unit)
    glucose_pattern = Chem.MolFromSmarts('C1OC(O)C(O)C(O)C(O)C1O')
    glucose_matches = mol.GetSubstructMatches(glucose_pattern)
    num_glucose_units = len(glucose_matches)

    # Check for ester linkages (-C(=O)O-)
    ester_pattern = Chem.MolFromSmarts('[CX3](=O)[OX2H0][#6]')
    num_ester_bonds = len(mol.GetSubstructMatches(ester_pattern))

    # Check for glycosidic bonds (oxygen linking two carbons in rings)
    glycosidic_pattern = Chem.MolFromSmarts('[$([C;R])]-O-[$([C;R])]')
    num_glycosidic_bonds = len(mol.GetSubstructMatches(glycosidic_pattern))

    # Determine if molecule meets criteria for tannin
    if num_glucose_units > 0 or num_glycosidic_bonds > 0 or num_ester_bonds > 0:
        reason = "Molecule has multiple catechol/pyrogallol units and linkages characteristic of tannins"
        return True, reason
    else:
        if num_phenol_groups >= 10:
            reason = "Molecule is a high molecular weight polyphenol with sufficient phenolic groups"
            return True, reason
        else:
            return False, "Insufficient linkages or sugar moieties characteristic of tannins"