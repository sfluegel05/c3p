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
        bool: True if molecule is a gamma-lactone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Identify gamma-lactone pattern: 5-membered lactone ring (O=C1OC(C)C1)
    gamma_lactone_pattern = Chem.MolFromSmarts("O=C1OC[C@H]2[C@H](O)CC2O1")
    if mol.HasSubstructMatch(gamma_lactone_pattern):
        return True, "Contains a 5-membered gamma-lactone ring"
    else:
        return False, "No gamma-lactone structure detected"

# Example test case
smiles_list = [
    "C/C=C/C1=C(C)CCCC1(C)C",  # Random example
    "O1CC(O)C(C1=O)=C"         # Beta-Hydroxy-alpha-methylene-gamma-butyllactone
]

for smiles in smiles_list:
    is_g_lactone, reason = is_gamma_lactone(smiles)
    print(f"SMILES: {smiles}, Gamma-Lactone: {is_g_lactone}, Reason: {reason}")