"""
Classifies: CHEBI:37581 gamma-lactone
"""
from rdkit import Chem

def is_gamma_lactone(smiles: str):
    """
    Determines if a molecule is a gamma-lactone based on its SMILES string.
    Gamma-lactones are characterized by a 5-membered lactone ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a gamma-lactone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Generalized pattern for a 5-membered gamma-lactone ring
    # Here [O] and [C] represent the basic elements of the ester linkage and ring carbons.
    # The bond to oxygen in the ring is specified as single bond, and the carbonyl is double bonded to oxygen.
    gamma_lactone_pattern = Chem.MolFromSmarts("C1C(=O)OC(C)C1")  # This represents a 5-membered gamma-lactone ring
    
    # Check if the molecule matches this pattern
    if mol.HasSubstructMatch(gamma_lactone_pattern):
        return True, "Contains a 5-membered gamma-lactone ring"
    else:
        return False, "No gamma-lactone structure detected"

# Example test cases for verification
smiles_list = [
    "O=C1OCCC1",  # A basic gamma-lactone structure
    "O=C1OC2(C)C(O)OCC12",  # Beta-Hydroxy-alpha-methylene-gamma-butyllactone
]

for smiles in smiles_list:
    is_g_lactone, reason = is_gamma_lactone(smiles)
    print(f"SMILES: {smiles}, Gamma-Lactone: {is_g_lactone}, Reason: {reason}")