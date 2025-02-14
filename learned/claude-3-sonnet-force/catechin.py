"""
Classifies: CHEBI:23053 catechin
"""
"""
Classifies: CHEBI:27776 catechin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_catechin(smiles: str):
    """
    Determines if a molecule is a catechin based on its SMILES string.
    A catechin is a member of the class of hydroxyflavan that has a flavan-3-ol skeleton and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a catechin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for flavonoid backbone pattern
    flavonoid_pattern = Chem.MolFromSmarts("[c&!r]1[C@H]2[C@@H]([C@H]([C@@H]([C@@]2(c2ccc(O)cc2)c2ccc(O)cc2)O)O)Oc3ccccc13")
    if not mol.HasSubstructMatch(flavonoid_pattern):
        flavonoid_pattern_alt = Chem.MolFromSmarts("[c&!r]1[C@@H]2[C@H]([C@H]([C@@H]([C@]2(c2ccc(O)cc2)c2ccc(O)cc2)O)O)Oc3ccccc13")
        if not mol.HasSubstructMatch(flavonoid_pattern_alt):
            return False, "No flavonoid backbone found"

    # Check for flavan-3-ol skeleton pattern
    flavan_pattern = Chem.MolFromSmarts("[C@@]1(c2ccccc2)Oc2ccccc2[C@@H]1[C@H](O)[C@H](O)[C@H](O)c1ccccc1")
    flavan_pattern_alt = Chem.MolFromSmarts("[C@]1(c2ccccc2)Oc2ccccc2[C@H]1[C@@H](O)[C@@H](O)[C@@H](O)c1ccccc1")
    if not mol.HasSubstructMatch(flavan_pattern) and not mol.HasSubstructMatch(flavan_pattern_alt):
        return False, "No flavan-3-ol skeleton found"

    # Check for common catechin substituents and modifications
    galloyl_pattern = Chem.MolFromSmarts("O=C(O)c1cc(O)c(O)c(O)c1")
    galloyl_matches = mol.GetSubstructMatches(galloyl_pattern)

    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    methoxy_pattern = Chem.MolFromSmarts("COc")
    methoxy_matches = mol.GetSubstructMatches(methoxy_pattern)

    hydroxyl_pattern = Chem.MolFromSmarts("O[H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)

    # Check for heterocyclic rings (pyran or chromene)
    pyran_pattern = Chem.MolFromSmarts("O1C=CC=CC1")
    chromene_pattern = Chem.MolFromSmarts("c1c2Oc3ccccc3Oc2ccc1")
    pyran_matches = mol.GetSubstructMatches(pyran_pattern)
    chromene_matches = mol.GetSubstructMatches(chromene_pattern)

    if not pyran_matches and not chromene_matches:
        return False, "No heterocyclic ring (pyran or chromene) found"

    # Check for relative positions of substituents (if possible)
    # (Omitted for brevity)

    # Count rings and aromatic rings
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    n_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if n_rings < 3 or n_rings > 6:
        return False, f"Number of rings ({n_rings}) out of expected range"
    if n_aromatic_rings < 2 or n_aromatic_rings > 4:
        return False, f"Number of aromatic rings ({n_aromatic_rings}) out of expected range"

    # Count hydroxyls (relaxed constraint)
    n_hydroxyls = len(hydroxyl_matches)
    if n_hydroxyls < 2:
        return False, f"Too few hydroxyl groups ({n_hydroxyls})"

    # Check molecular formula pattern
    mol_formula = rdMolDescriptors.CalcMolFormula(mol)
    catechin_formula_pattern = r"C\d{12,24}H\d{14,30}O\d{4,12}"
    if not re.match(catechin_formula_pattern, mol_formula):
        return False, f"Molecular formula ({mol_formula}) does not match expected pattern"

    # If all checks pass, classify as catechin
    return True, "Molecule satisfies structural and chemical features of catechins"