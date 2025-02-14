"""
Classifies: CHEBI:26493 quinic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is a quinic acid based on its SMILES string.
    A quinic acid is a cyclitol carboxylic acid, having a cyclohexane ring with multiple hydroxyl groups and a carboxyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Core structure: cyclohexane with COOH. At least 3 oxygens attached to ring
    # Modified to find the basic structure with a COOH and at least 3 oxygen atoms attached
    quinic_acid_pattern = Chem.MolFromSmarts("[CX4]1([CX4]([OX2])([CX4])([CX4])([CX4])([CX4]1))C(=O)[OX1H0]")
    if mol.HasSubstructMatch(quinic_acid_pattern):
        # Check that the ring has at least 3 oxygen atoms
        ring_atoms = [atom for atom in mol.GetAtoms() if atom.IsInRing() and atom.GetAtomicNum() == 6]
        oxygens_on_ring = 0
        for ring_carbon in ring_atoms:
            for neighbor in ring_carbon.GetNeighbors():
                if neighbor.GetAtomicNum() == 8:
                    oxygens_on_ring+=1
        if oxygens_on_ring >= 3:
            return True, "Matches quinic acid pattern (cyclohexane ring with at least 3 oxygen atoms and COOH)"

    return False, "Does not match quinic acid pattern"