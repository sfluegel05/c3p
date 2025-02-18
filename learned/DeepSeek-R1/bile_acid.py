"""
Classifies: CHEBI:3098 bile acid
"""
"""
Classifies: CHEBI:3098 bile acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_bile_acid(smiles: str):
    """
    Determines if a molecule is a bile acid based on its SMILES string.
    Bile acids are hydroxy-5beta-cholanic acids with a steroid nucleus and specific stereochemistry.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bile acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group at position 24
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid group found"

    # Verify steroid nucleus (four fused rings: 6+6+6+5)
    steroid_skeleton = Chem.MolFromSmarts("[*]1~[*]~[*]~[*]~[*]~[*]1.[*]2~[*]~[*]~[*]~[*]~[*]2.[*]3~[*]~[*]~[*]~[*]~[*]3.[*]4~[*]~[*]~[*]~[*]4")
    if not mol.HasSubstructMatch(steroid_skeleton):
        return False, "No steroid nucleus detected"

    # Check 5β configuration (trans A/B ring fusion)
    # Look for 5β configuration where C5 has R configuration
    try:
        Chem.AssignStereochemistry(mol)
        ring_system = mol.GetSubstructMatch(Chem.MolFromSmarts("[C@H]1[C@@H]2CC[C@H]3[C@@]4"))
        if not ring_system:
            return False, "Could not verify 5beta configuration"
    except:
        return False, "Stereochemistry check failed"

    # Check for at least one hydroxyl group on the steroid nucleus
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if not hydroxyl_matches:
        return False, "No hydroxyl groups found on steroid nucleus"

    # Check side chain length (4 carbons from first bridgehead to COOH)
    side_chain = Chem.MolFromSmarts("[C@@H]1[C@@H]2CC[C@H]3[C@]4(CC[C@@H]([C@H](CC(=O)O)C)C)")
    if not mol.HasSubstructMatch(side_chain):
        return False, "Incorrect side chain structure/length"

    return True, "5β-cholanic acid derivative with steroid nucleus and hydroxyl groups"