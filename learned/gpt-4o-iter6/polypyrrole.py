"""
Classifies: CHEBI:38077 polypyrrole
"""
from rdkit import Chem

def is_polypyrrole(smiles: str):
    """
    Determines if a molecule is a polypyrrole based on its SMILES string.
    A polypyrrole is a compound composed of two or more pyrrole units or pyrrole-like structures.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        tuple: A tuple containing a boolean indicating if the molecule is a polypyrrole,
               and a reason for the classification.
    """
    # Parse SMILES string to RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a comprehensive list of SMARTS for pyrrole and polypyrrole structures
    pyrrole_patterns = [
        Chem.MolFromSmarts("c1c[nH]c1"),  # Basic Pyrrole
        Chem.MolFromSmarts("n1cccc1"),    # Pyrrole
        Chem.MolFromSmarts("n1ccccc1"),   # Larger systems incorporating pyrrole
        Chem.MolFromSmarts("c1c[nH]c2c1nccc2"),  # Fused pyrrole systems
        Chem.MolFromSmarts("n2c1nc(nc1[nH]2)"),  # Azole-like structures
        Chem.MolFromSmarts("c1ccc2c(c1)nccc2n"), # Expanded porphyrin-like structure
    ]
    
    # Count all occurrences of pyrrole-like structures
    num_pyrrole_structures = 0
    for pattern in pyrrole_patterns:
        if pattern is None:
            continue

        # Find matches for each pattern and accumulate the count
        matches = mol.GetSubstructMatches(pattern)
        num_pyrrole_structures += len(matches)

    # Check for polypyrrole condition: at least two pyrrole or pyrrole-like structures
    if num_pyrrole_structures >= 2:
        return True, f"Contains {num_pyrrole_structures} pyrrole-like structures"
    else:
        return False, "Contains less than two pyrrole-like structures"

# Example usage:
smiles_example = "CCC1=C(C)C(=O)NC1Cc1[nH]c(Cc2[nH]c(CC3NC(=O)C(CC)=C3C)c(C)c2CCC(O)=O)c(CCC(O)=O)c1C"
print(is_polypyrrole(smiles_example))