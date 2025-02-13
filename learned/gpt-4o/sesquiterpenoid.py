"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a sesquiterpenoid based on its SMILES string.
    A sesquiterpenoid is any terpenoid derived from a sesquiterpene, often modified by rearrangement 
    or methyl group removal, having a C15 carbon skeleton.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a sesquiterpenoid, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of carbon skeleton typically related to sesquiterpenoids
    sesquiterpene_pattern = Chem.MolFromSmarts('C1C(C)CCC1')
    if not mol.HasSubstructMatch(sesquiterpene_pattern):
        return False, "Does not contain a substructure typical of sesquiterpenoid backbone"

    # Check for presence of oxygen indicating functionalization
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 1:
        return False, "Does not contain oxygen, unlikely to be a sesquiterpenoid"

    # Check for additional complex sesquiterpenoid structures or functional groups
    complex_sesquiterpenoid_patterns = [
        Chem.MolFromSmarts('C=1C[CH]C2=C1[CH][CH2]C=CC2')  # Example complex pattern smarts
    ]

    for pattern in complex_sesquiterpenoid_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains complex sesquiterpenoid pattern"

    # If passes backbone and modification checks, it's likely a sesquiterpenoid
    return True, "Contains characteristics typical of a sesquiterpenoid"

# Example
print(is_sesquiterpenoid('O=C1C=C2C=CC(=O)[C@@]([C@]2(C)C[C@]1(O)C(=C)COC(=O)C)(O)C'))  # Expected True classification