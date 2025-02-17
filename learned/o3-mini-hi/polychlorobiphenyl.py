"""
Classifies: CHEBI:53156 polychlorobiphenyl
"""
"""
Classifies: Polychlorobiphenyl
Definition: A biphenyl compound containing between 2 and 10 chlorine atoms attached to the two benzene rings.
"""
from rdkit import Chem

def is_polychlorobiphenyl(smiles: str):
    """
    Determines if a molecule is a polychlorobiphenyl based on its SMILES string.
    
    A polychlorobiphenyl must have:
         - A biphenyl scaffold: two benzene rings (aromatic six-membered rings) connected by a single bond.
         - Between 2 and 10 chlorine substituents.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is a polychlorobiphenyl, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for the biphenyl scaffold.
    # This pattern looks for two connected benzene rings.
    biphenyl_pattern = Chem.MolFromSmarts("c1ccc(cc1)-c2ccccc2")
    if not mol.HasSubstructMatch(biphenyl_pattern):
        return False, "Biphenyl scaffold not detected"
    
    # Count the number of chlorine (Cl) atoms in the molecule.
    chlorine_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 17)
    if chlorine_count < 2:
        return False, f"Found {chlorine_count} chlorine atoms; at least 2 are required"
    if chlorine_count > 10:
        return False, f"Found {chlorine_count} chlorine atoms; no more than 10 are allowed"
    
    return True, f"Contains biphenyl scaffold with {chlorine_count} chlorine atoms"