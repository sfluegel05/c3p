"""
Classifies: CHEBI:10615 omega-hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdChemReactions

def is_omega_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an omega-hydroxy fatty acid based on its SMILES string.
    An omega-hydroxy fatty acid has a carboxyl group at position 1 and a hydroxyl at position n (omega).
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is an omega-hydroxy fatty acid, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern for carboxylic acid
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not carboxyl_matches:
        return False, "No carboxyl group present"
    
    # Pattern for terminal hydroxyl (omega position)
    hydroxyl_pattern = Chem.MolFromSmarts("O")
    hydroxyl_matches = [match for match in mol.GetSubstructMatches(hydroxyl_pattern) if len(match) == 1 and mol.GetAtomWithIdx(match[0]).GetDegree() == 1]
    if not hydroxyl_matches:
        return False, "No terminal hydroxyl group found at omega position"
    
    # Count carbon chain length
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 6:
        return False, f"Carbon chain too short has only {carbon_count} carbons"

    # Verify distinct ends: one carboxyl end and one omega -OH end
    if len(carboxyl_matches) != 1 or len(hydroxyl_matches) != 1:
        return False, "Should contain one carboxyl group and one omega-hydroxyl group"
    
    return True, "Contains carboxyl group and omega-hydroxyl group as per omega-hydroxy fatty acid definition."