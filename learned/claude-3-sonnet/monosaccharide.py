"""
Classifies: CHEBI:35381 monosaccharide
"""
"""
Classifies: CHEBI:16646 monosaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monosaccharide(smiles: str):
    """
    Determines if a molecule is a monosaccharide based on its SMILES string.
    Parent monosaccharides are polyhydroxy aldehydes H[CH(OH)]nC(=O)H or polyhydroxy ketones
    H-[CHOH]n-C(=O)[CHOH]m-H with three or more carbon atoms. The generic term 'monosaccharide'
    (as opposed to oligosaccharide or polysaccharide) denotes a single unit, without glycosidic
    connection to other such units. It includes aldoses, dialdoses, aldoketoses, ketoses and
    diketoses, as well as deoxy sugars, provided that the parent compound has a (potential) carbonyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monosaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carbonyl group (aldehyde or ketone)
    aldehyde_pattern = Chem.MolFromSmarts("C=O")
    ketone_pattern = Chem.MolFromSmarts("C(=O)C")
    has_carbonyl = mol.HasSubstructMatch(aldehyde_pattern) or mol.HasSubstructMatch(ketone_pattern)
    if not has_carbonyl:
        return False, "No carbonyl group (aldehyde or ketone) found"

    # Check for hydroxyl groups (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX1H]")
    n_hydroxyls = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if n_hydroxyls < 3:
        return False, "Fewer than 3 hydroxyl groups found"

    # Check for monosaccharide skeletons (aldoses, ketoses, furanoses, pyranoses)
    aldose_pattern = Chem.MolFromSmarts("C(C(O)O)O")
    ketose_pattern = Chem.MolFromSmarts("C(C(O)O)(C(O)O)O")
    furanose_pattern = Chem.MolFromSmarts("C1OC(O)C(O)C1O")
    pyranose_pattern = Chem.MolFromSmarts("C1OC(O)C(O)C(O)C(O)C1O")
    has_monosaccharide_skeleton = (
        mol.HasSubstructMatch(aldose_pattern) or
        mol.HasSubstructMatch(ketose_pattern) or
        mol.HasSubstructMatch(furanose_pattern) or
        mol.HasSubstructMatch(pyranose_pattern)
    )
    if not has_monosaccharide_skeleton:
        return False, "No monosaccharide skeleton found"

    # Check for glycosidic bonds (to exclude oligosaccharides and polysaccharides)
    glycosidic_pattern = Chem.MolFromSmarts("O[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)O")
    has_glycosidic_bond = mol.HasSubstructMatch(glycosidic_pattern)
    if has_glycosidic_bond:
        return False, "Glycosidic bond found, not a single monosaccharide unit"

    # Check molecular weight and size constraints
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    n_atoms = mol.GetNumAtoms()
    if mol_wt < 100 or mol_wt > 400 or n_atoms < 5 or n_atoms > 20:
        return False, "Molecular weight or size outside typical range for monosaccharides"

    return True, "Molecule matches the structure of a monosaccharide"