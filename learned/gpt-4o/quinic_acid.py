"""
Classifies: CHEBI:26493 quinic acid
"""
from rdkit import Chem

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is a quinic acid or its derivative based on its SMILES string.
    A quinic acid derivative typically consists of a cyclohexane core with hydroxy and ester groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinic acid or relevant derivative, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for a cyclohexane ring
    ring_info = mol.GetRingInfo()
    if not any(len(ring) == 6 for ring in ring_info.AtomRings()):
        return False, "No cyclohexane ring found"

    # Check for hydroxy and ester groups on the cyclohexane, allowing flexible stereochemistry
    hydroxy_pattern = Chem.MolFromSmarts("C(O)")
    ester_pattern = Chem.MolFromSmarts("OC=O")
    
    # Check complex quinic features within the ring
    pattern = Chem.MolFromSmarts("[C;R1]1([OH])[C;R1][C;R1](O)[C;R1]([OH])[C;R1](OC=O)[C;R1]1C(=O)O")
    
    if not mol.HasSubstructMatch(pattern):
        # Check if it matches a simpler pattern where stereochemistry might not be resolved
        rearranged_pattern = Chem.MolFromSmarts("[C;R1]1([OH])[C;R1][C;R1](O)[C;R1][C;R1](O)[C;R1]1C(=O)O")
        if not mol.HasSubstructMatch(rearranged_pattern):
            return False, "Does not match quinic acid core patterns"
    
    if any(map(mol.HasSubstructMatch, [hydroxy_pattern, ester_pattern])):
        return True, "Contains features consistent with esterified quinic acid"

    return False, "Features do not match quinic acid derivatives"

# Example use case
smiles_example = "O[C@H]1C[C@@](O)(C[C@H](O)[C@H]1O)C(O)=O"
print(is_quinic_acid(smiles_example))