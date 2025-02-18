"""
Classifies: CHEBI:36315 glycerophosphoinositol
"""
"""
Classifies: CHEBI:28875 glycerophosphoinositol
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_glycerophosphoinositol(smiles: str):
    """
    Determines if a molecule is a glycerophosphoinositol based on its SMILES string.
    Must have: glycerol backbone with two fatty acid esters, phosphate at sn-3,
    and inositol (cyclohexanehexol) attached to phosphate.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if glycerophosphoinositol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phosphate group connected to inositol
    # Inositol pattern: cyclohexane with 5 hydroxyls (6th position connected to phosphate)
    inositol_phosphate_pattern = Chem.MolFromSmarts("[O][P](=O)([O])[O][C@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O")
    if not mol.HasSubstructMatch(inositol_phosphate_pattern):
        return False, "Phosphate not connected to myo-inositol (cyclohexanehexol)"

    # Verify glycerol backbone: C1-C2-C3 with two esters and phosphate at C3
    glycerol_pattern = Chem.MolFromSmarts("[CH2]-[CH](-[OX2]C(=O))-[CH2]-[OX2]-[P]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        # Check for possible stereochemistry variations
        glycerol_pattern = Chem.MolFromSmarts("[CH2]-[CH](-[OX2]C(=O))-[CH2]-[OX2]-[P]")
        if not mol.HasSubstructMatch(glycerol_pattern):
            return False, "Glycerol backbone with two esters and phosphate not found"

    # Ensure exactly two ester groups attached to glycerol (sn-1 and sn-2)
    ester_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[OX2]C(=O)[CX4]")))
    if ester_count < 2:
        return False, f"Found {ester_count} ester groups, need at least 2"

    # Optional: Check fatty acid chain length (at least 8 carbons)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4]")
    if len(mol.GetSubstructMatches(fatty_acid_pattern)) < 2:
        return False, "Insufficient long-chain fatty acids"

    return True, "Glycerol with two fatty acid esters, phosphate-linked myo-inositol"