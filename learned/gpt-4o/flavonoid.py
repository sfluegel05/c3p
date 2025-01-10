"""
Classifies: CHEBI:47916 flavonoid
"""
from rdkit import Chem

def is_flavonoid(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    A flavonoid is a member of the superclass with a 1-benzopyran understructure and variations 
    of aryl substituent at position 2 with possible multiple substitutions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern capturing 1-benzopyran and variations with flexibility for modification
    flavonoid_smarts = "[O]c1cc(-c2ccccc2)cc2c1cc(O)cc2"
    flavonoid_pattern = Chem.MolFromSmarts(flavonoid_smarts)

    # Check for flavonoid structure
    if not mol.HasSubstructMatch(flavonoid_pattern):
        return False, "No flavonoid structure detected"

    return True, "Molecule contains the 1-benzopyran structure with aryl group typical of flavonoids"

# Test with sample SMILES demonstrating flavonoid structure
sample_smiles = [
    "O1C2=C(C(OC)=C3C(OC=C3)=C2)C(=O)C=C1C4=CC=CC=C4",  # Pinnatin
    "O1C2=C(CC=C(C)C)C(O)=CC(O)=C2C(=O)C=C1C3=CC(OC)=C(O)C(OC)=C3",  # Baohuosu
]

results = [(smiles, is_flavonoid(smiles)) for smiles in sample_smiles]
for smiles, (result, reason) in results:
    print(f"SMILES: {smiles}\nResult: {result}\nReason: {reason}\n")