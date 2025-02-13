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

    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    # Check for typical sesquiterpenoid carbon count (C15)
    if c_count < 13 or c_count > 15:
        return False, f"Carbon count is {c_count}, not typical for sesquiterpenoid (C13 to C15 after modifications)"

    # General terpenoid characteristic: presence of oxygen
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 1:
        return False, "Does not contain oxygen, unlikely to be a sesquiterpenoid"

    # Optional: identify rearrangements and potential skeletal modifications
    # This could involve specific substructure searches typical of known sesquiterpenoid scaffolds

    # If all checks pass, classify as sesquiterpenoid
    return True, "Contains characteristics typical of a sesquiterpenoid (C15 related skeleton with functional groups such as oxygen)"

# Example
print(is_sesquiterpenoid('O=C1C=C2C=CC(=O)[C@@]([C@]2(C)C[C@]1(O)C(=C)COC(=O)C)(O)C'))  # Output will be a tuple indicating classification