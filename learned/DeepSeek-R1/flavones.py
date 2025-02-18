"""
Classifies: CHEBI:24043 flavones
"""
"""
Classifies: CHEBI:72025 flavone
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_flavones(smiles: str):
    """
    Determines if a molecule is a flavone based on its SMILES string.
    Flavones have a 2-aryl-1-benzopyran-4-one skeleton and substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for the 2-arylchromen-4-one skeleton
    # Chromen-4-one core (benzopyran-4-one) with aryl group at position 2
    flavone_pattern = Chem.MolFromSmarts("c1c2c(oc(=O)cc2)ccc1-[a]")
    if not mol.HasSubstructMatch(flavone_pattern):
        return False, "No 2-arylchromen-4-one skeleton found"

    # Additional check to ensure the aryl group is a benzene ring (optional)
    # This can be adjusted based on specific requirements
    # Get the matching atoms
    matches = mol.GetSubstructMatches(flavone_pattern)
    for match in matches:
        # The last atom in the match is the aryl group's connecting atom
        aryl_connector = match[-1]
        # Check if it's part of an aromatic ring (benzene or substituted)
        neighbor = mol.GetAtomWithIdx(aryl_connector)
        if neighbor.GetIsAromatic():
            return True, "Contains 2-arylchromen-4-one skeleton"

    return False, "Aryl group not aromatic or not properly attached"