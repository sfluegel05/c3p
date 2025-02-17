"""
Classifies: CHEBI:26377 pterocarpans
"""
"""
Classifies: Pterocarpans – Members of the class of benzofurochromene with a 6a,11a-dihydro-6H-[1]benzofuro[3,2-c]chromene skeleton (and substituted derivatives).
They are 3,4-dihydro derivatives of coumestans and generally resemble isoflavanoids.
"""

from rdkit import Chem

def is_pterocarpans(smiles: str):
    """
    Determines if a molecule is a pterocarpan based on its SMILES string.
    
    The classification is based on the presence of the benzofurochromene core.
    A simplified core is represented by a dihydrofuran ring fused to two benzene rings.
    The SMARTS pattern used here is "C1COc2ccccc2Oc3ccccc13", which captures this tricyclic system.
    Note: This pattern is a simplified representation and will capture many substituted derivatives.
    
    Args:
        smiles (str): SMILES representation of the molecule.

    Returns:
        bool: True if the molecule contains the pterocarpan core, False otherwise.
        str: Explanation of the classification result.
    """
    # Parse the SMILES string into an RDKit Mol object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern corresponding to the pterocarpan core.
    # The pattern "C1COc2ccccc2Oc3ccccc13" represents a dihydrofuran (C1CO)
    # fused on one side to a benzene ring (c2ccccc2) and on the other to another benzene ring (c3ccccc1 with closure).
    pattern = Chem.MolFromSmarts("C1COc2ccccc2Oc3ccccc13")
    if pattern is None:
        return False, "Internal error: pattern could not be parsed"
    
    # Check if the pterocarpan core is present in the molecule.
    if not mol.HasSubstructMatch(pattern):
        return False, "Pterocarpan core (benzofurochromene skeleton) not found"
    
    return True, "Molecule contains the pterocarpan core scaffold"