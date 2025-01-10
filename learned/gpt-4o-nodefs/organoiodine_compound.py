"""
Classifies: CHEBI:37142 organoiodine compound
"""
from rdkit import Chem

def is_organoiodine_compound(smiles: str):
    """
    Determines if a molecule is an organoiodine compound based on its SMILES string.
    An organoiodine compound is characterized by broader structural features typically
    involving iodine covalently bonded to carbon, often in specific arrangements.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organoiodine compound, False otherwise
        str: Reason for classification or misclassification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for iodine atoms within the molecule bonded directly to carbon
    iodine_carbon_bond_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 53:  # Atomic number of iodine is 53
            # Check if this iodine is bonded directly to a carbon atom
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:  # Atomic number of carbon is 6
                    iodine_carbon_bond_found = True
                    break
            if iodine_carbon_bond_found:
                break

    # If iodine-carbon bond found, additional structure checks can be added here,
    # like verifying the presence of common organoiodine substructures.
    # Placeholder for more elaborate pattern matching if required later:
    # Example:
    # organoiodine_pattern = Chem.MolFromSmarts('<specific_SMARTS_here>')
    # if not mol.HasSubstructMatch(organoiodine_pattern):
    #     return False, f"Pattern common to organoiodine compounds not found."

    # For demonstration, consider iodine-carbon bond presence as a marker
    if iodine_carbon_bond_found:
        return True, "Iodine is bonded directly to carbon, indicating an organoiodine compound."

    return False, "No iodine-carbon bond found; not an organoiodine compound."