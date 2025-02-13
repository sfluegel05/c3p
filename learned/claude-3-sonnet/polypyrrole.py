"""
Classifies: CHEBI:38077 polypyrrole
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_polypyrrole(smiles: str):
    """
    Determines if a molecule is a polypyrrole based on its SMILES string.
    A polypyrrole is a compound composed of two or more pyrrole units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polypyrrole, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for a connected pyrrole unit
    pyrrole_pattern = Chem.MolFromSmarts("[nH]1ccc[nH]1")
    
    # Find all matches of the pyrrole pattern in the molecule
    pyrrole_matches = mol.GetSubstructMatches(pyrrole_pattern)
    
    # Keep track of the total number of pyrrole units
    total_pyrrole_units = 0
    
    # Iterate over matches to identify connected pyrrole fragments
    for match in pyrrole_matches:
        pyrrole_fragment = []
        for atom_idx in match:
            pyrrole_fragment.append(atom_idx)
            # Check for connected pyrrole units
            connected_pyrroles = mol.GetAtomWithIdx(atom_idx).GetNeighbors()
            for neighbor in connected_pyrroles:
                if neighbor.GetIdx() not in pyrrole_fragment and mol.GetAtomWithIdx(neighbor.GetIdx()).GetIsAromaticHeteroatom():
                    pyrrole_fragment.append(neighbor.GetIdx())
        
        # Count the number of pyrrole units in the fragment
        pyrrole_units_in_fragment = len(pyrrole_fragment) // 5
        total_pyrrole_units += pyrrole_units_in_fragment
    
    # Classify as a polypyrrole if at least two pyrrole units are found
    if total_pyrrole_units >= 2:
        return True, f"Contains {total_pyrrole_units} connected pyrrole units"
    else:
        return False, "Less than two connected pyrrole units found"