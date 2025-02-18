"""
Classifies: CHEBI:26493 quinic acid
"""
"""
Classifies: CHEBI:17544 quinic acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is a quinic acid derivative based on its SMILES string.
    A quinic acid is a cyclitol carboxylic acid with a cyclohexane backbone, multiple hydroxyl groups,
    and a carboxylic acid group. Derivatives include esters of these hydroxyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinic acid derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-COOH)
    carboxylic_acid = Chem.MolFromSmarts('[CX3](=O)[OX2H1]')
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid group present"

    # Find all 6-membered carbon rings
    ring_info = mol.GetRingInfo()
    six_membered_carbon_rings = [
        ring for ring in ring_info.AtomRings()
        if len(ring) == 6 and all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring)
    ]
    if not six_membered_carbon_rings:
        return False, "No six-membered carbon ring found"

    # Check if carboxylic acid is attached to the ring
    # SMARTS for ring carbon connected to COOH: [C;R]C(=O)O
    ring_carboxylic = Chem.MolFromSmarts('[C;R]C(=O)O')
    if not mol.HasSubstructMatch(ring_carboxylic):
        return False, "Carboxylic acid not attached to cyclohexane ring"

    # Count oxygen atoms attached to ring carbons (including those in esters and hydroxyls)
    total_ring_oxygens = 0
    for ring in six_membered_carbon_rings:
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            total_ring_oxygens += sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 8)

    # Quinic acid core has at least 4 oxygens (COOH contributes 2, others contribute at least 2)
    if total_ring_oxygens < 4:
        return False, f"Insufficient oxygen substituents on ring: {total_ring_oxygens}"

    return True, "Cyclohexane carboxylic acid with sufficient oxygen substituents"