"""
Classifies: CHEBI:61703 nonclassic icosanoid
"""
"""
Classifies: nonclassic icosanoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nonclassic_icosanoid(smiles: str):
    """
    Determines if a molecule is a nonclassic icosanoid based on its SMILES string.
    Nonclassic icosanoids are biologically active C20 fatty acid derivatives with oxygenation, excluding prostanoids and leukotrienes.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check total carbons = 20 (derived from C20 fatty acid)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 20:
        return False, f"Contains {c_count} carbons (requires 20)"

    # Check for carboxylic acid or ester group (common in icosanoids)
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")  # Carboxylic acid
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2][#6]")  # Ester
    if not mol.HasSubstructMatch(acid_pattern) and not mol.HasSubstructMatch(ester_pattern):
        return False, "No carboxylic acid or ester group"

    # Check oxygenation: at least one hydroxyl, epoxy, or ketone (excluding acid/ester)
    hydroxyl = mol.HasSubstructMatch(Chem.MolFromSmarts("[OH]"))
    epoxy = mol.HasSubstructMatch(Chem.MolFromSmarts("C1OC1"))
    ketone = mol.HasSubstructMatch(Chem.MolFromSmarts("[CX3]=[OX1]"))  # Exclude acid/ester (already checked)
    if not (hydroxyl or epoxy or ketone):
        return False, "No hydroxyl, epoxy, or ketone groups"

    # Exclude prostanoids (5-membered carbon ring)
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) == 5 and all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
            return False, "Contains 5-membered carbon ring (prostanoid)"

    # Exclude leukotrienes (conjugated trienes: three consecutive double bonds)
    triene_pattern = Chem.MolFromSmarts("C=C-C=C-C=C")
    if mol.HasSubstructMatch(triene_pattern):
        return False, "Contains conjugated triene (leukotriene)"

    return True, "C20 fatty acid derivative with oxygenation, not prostanoid/leukotriene"