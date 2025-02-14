"""
Classifies: CHEBI:51129 nitrohydrocarbon
"""
from rdkit import Chem

def is_nitrohydrocarbon(smiles: str):
    """
    Determines if a molecule is a nitrohydrocarbon based on its SMILES string.
    A nitrohydrocarbon is a hydrocarbon in which one or more of the hydrogens has been replaced by nitro groups.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nitrohydrocarbon, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the nitro group pattern
    nitro_group = Chem.MolFromSmarts("[NX3](=O)[O-]")
    
    # Check for nitro groups attached to carbon
    found_nitro_on_carbon = False
    for match in mol.GetSubstructMatches(nitro_group):
        # Check if the nitro N is connected to a carbon atom
        nitro_nitrogen = mol.GetAtomWithIdx(match[0])
        if any(neighbor.GetAtomicNum() == 6 for neighbor in nitro_nitrogen.GetNeighbors()):
            found_nitro_on_carbon = True
            break
    
    if not found_nitro_on_carbon:
        return False, "No nitro group replacing hydrogen on a carbon atom found"
    
    # Ensure the molecule is primarily a hydrocarbon backbone
    permissible_atoms = {6, 1, 7, 8}  # C, H, N (in nitro), O (in nitro)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in permissible_atoms:
            return False, f"Contains non-hydrocarbon atoms: {atom.GetSymbol()}"
   
    return True, "Contains hydrocarbon structure with one or more nitro groups replacing hydrogen"