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
    Tannins are polyphenolic compounds with multiple galloyl groups, ellagitannin features,
    or condensed flavan-3-ol structures.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tannin, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Improved galloyl group detection (3,4,5-trihydroxybenzoate ester)
    galloyl_pattern = Chem.MolFromSmarts('[OX2]C(=O)c1c(O)c(O)c(O)cc1')  # Corrected positions
    galloyl_count = len(mol.GetSubstructMatches(galloyl_pattern))
    if galloyl_count >= 1:
        # Check for multiple galloyl esters or presence with other phenolic groups
        ester_pattern = Chem.MolFromSmarts('[OX2]C(=O)')  # General ester count
        ester_count = len(mol.GetSubstructMatches(ester_pattern))
        if ester_count >= 2 and galloyl_count >= 1:
            return True, f"Contains {galloyl_count} galloyl groups and {ester_count} esters"
        if galloyl_count >= 2:
            return True, f"Contains {galloyl_count} galloyl groups"

    # Check for hexahydroxydiphenoyl (HHDP) groups (characteristic of ellagitannins)
    hhdp_pattern = Chem.MolFromSmarts('c1c(O)c(O)c(-c2c(O)c(O)c(O)cc2)c(O)c1')
    if mol.HasSubstructMatch(hhdp_pattern):
        return True, "Contains hexahydroxydiphenoyl (HHDP) group"

    # Detect flavan-3-ol oligomers (condensed tannins)
    # Basic flavan-3-ol pattern with multiple hydroxyls
    flavan_pattern = Chem.MolFromSmarts('[C@H]1Oc2cc(O)cc(O)c2C[C@H](O)C1')
    flavan_matches = len(mol.GetSubstructMatches(flavan_pattern))
    if flavan_matches >= 2:
        return True, f"Contains {flavan_matches} flavan-3-ol subunits"

    # Enhanced phenolic group detection with positional requirements
    pyrogallol = Chem.MolFromSmarts('c1c(O)c(O)c(O)cc1')  # Three adjacent OH
    catechol = Chem.MolFromSmarts('c1c(O)c(O)ccc1')       # Two adjacent OH
    
    phenolic_count = len(mol.GetSubstructMatches(pyrogallol)) * 3 + len(mol.GetSubstructMatches(catechol)) * 2
    total_oh = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() >= 1)
    
    # Check for high hydroxylation and molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if (phenolic_count >= 6 or total_oh >= 8) and mol_wt > 500:
        return True, f"Highly hydroxylated ({total_oh} OH) with MW {mol_wt:.1f}"

    # Check for glucose core with multiple substituents (common in hydrolyzable tannins)
    glucose_core = Chem.MolFromSmarts('[C@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O')
    if mol.HasSubstructMatch(glucose_core):
        # Look for at least 3 ester groups attached to glucose
        ester_attachments = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[OX2]C(=O)-[#6]@;-[OX2][C@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O')))
        if ester_attachments >= 3:
            return True, f"Glucose core with {ester_attachments} ester groups"

    return False, "Does not meet tannin criteria"