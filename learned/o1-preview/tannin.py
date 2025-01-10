"""
Classifies: CHEBI:26848 tannin
"""
"""
Classifies: CHEBI:27027 tannin
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tannin(smiles: str):
    """
    Determines if a molecule is a tannin based on its SMILES string.
    Tannins are a group of astringent polyphenolic compounds, chiefly complex glucosides of catechol and pyrogallol.
    They include hydrolyzable tannins (gallotannins and ellagitannins) and condensed tannins (proanthocyanidins).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tannin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns

    # Gallic acid unit: a phenolic ring with hydroxyls at positions 3, 4, and 5
    gallic_acid_pattern = Chem.MolFromSmarts('c1c(O)cc(O)c(O)c1')

    # Ellagic acid unit: biphenyl structure with multiple hydroxyls and lactone rings
    ellagic_acid_pattern = Chem.MolFromSmarts('O=C1c2c(O)c(O)cc(O)c2C(=O)c2c1c(O)c(O)cc2O')

    # Flavan-3-ol unit (catechin/epicatechin): core unit of condensed tannins
    flavan3ol_pattern = Chem.MolFromSmarts('C1[C@@H](O)[C@@H](Oc2cc(O)cc(O)c12)c1cc(O)cc(O)c1')

    # Interflavan linkage C4-C8 (common in proanthocyanidins)
    interflavan_C4_C8_pattern = Chem.MolFromSmarts('[C@@H]1(O)[C@H](Oc2ccc(O)cc2)[C@H](O[C@@H]3[C@@H](O)[C@H](Oc4ccc(O)cc34)C1)c1ccc(O)cc1')

    # Ester bond (for gallotannins and ellagitannins)
    ester_pattern = Chem.MolFromSmarts('C(=O)O[C;!$(C=O)]')

    # Glycosidic linkage: sugar linkage (C-O-C between carbons)
    glycosidic_pattern = Chem.MolFromSmarts('[#6]-O-[#6]')

    # Count matches of key substructures
    num_gallic_acid = len(mol.GetSubstructMatches(gallic_acid_pattern))
    num_ellagic_acid = len(mol.GetSubstructMatches(ellagic_acid_pattern))
    num_flavan3ol = len(mol.GetSubstructMatches(flavan3ol_pattern))
    num_interflavan = len(mol.GetSubstructMatches(interflavan_C4_C8_pattern))
    num_esters = len(mol.GetSubstructMatches(ester_pattern))
    num_glycosidic = len(mol.GetSubstructMatches(glycosidic_pattern))

    # Total phenolic units (gallic acid units, ellagic acid units, flavan-3-ol units)
    total_phenolic_units = num_gallic_acid + num_ellagic_acid + num_flavan3ol

    # Check for criteria of tannins

    if total_phenolic_units < 2:
        return False, f"Found {total_phenolic_units} phenolic units, need at least 2"

    # Tannins often have multiple ester or glycosidic bonds
    if num_esters + num_glycosidic < 1:
        return False, f"Found {num_esters} ester bonds and {num_glycosidic} glycosidic bonds, need at least 1 combined"

    # High degree of polymerization (aromatic rings >= 4)
    num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if num_aromatic_rings < 3:
        return False, f"Found {num_aromatic_rings} aromatic rings, need at least 3 for tannins"

    # All criteria met
    return True, "Contains multiple phenolic units with ester or glycosidic bonds and sufficient aromatic rings characteristic of tannins"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:27027',
        'name': 'tannin',
        'definition': 'Any of a group of astringent polyphenolic vegetable principles or compounds, chiefly complex glucosides of catechol and pyrogallol.',
        'parents': ['CHEBI:26195']
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
        'f1_threshold': 0.8,
        'max_attempts': 5,
        'max_positive_instances': None,
        'max_positive_to_test': None,
        'max_negative_to_test': None,
        'max_positive_in_prompt': 50,
        'max_negative_in_prompt': 20,
        'max_instances_in_prompt': 100,
        'test_proportion': 0.1
    }
}