"""
Classifies: CHEBI:143084 organometalloidal compound
"""
from rdkit import Chem

def is_organometalloidal_compound(smiles: str):
    """
    Determines if a molecule is an organometalloidal compound based on the presence of metalloid-carbon bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple: (bool, str) where bool is True if organometalloidal, False otherwise;
               str provides the reason. Returns (None, None) for invalid SMILES.
    """
    # Define metalloids by atomic number (As, B, Si, Ge, Sb, Te, Po)
    metalloids = {5, 14, 32, 33, 51, 52, 84}
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None  # Invalid SMILES
    
    # Iterate through each atom in the molecule
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in metalloids:
            # Check all neighboring atoms for carbon (atomic number 6)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:
                    return True, f"Contains {atom.GetSymbol()}-C bond"
    
    # No metalloid-carbon bonds found
    return False, "No metalloid-carbon bonds detected"