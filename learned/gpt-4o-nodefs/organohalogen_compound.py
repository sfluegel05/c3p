"""
Classifies: CHEBI:17792 organohalogen compound
"""
"""
Classifies: Organohalogen compounds
"""
from rdkit import Chem

def is_organohalogen_compound(smiles: str):
    """
    Determines if a molecule is an organohalogen compound based on its SMILES string.
    Organohalogen compounds have carbon atoms bonded to halogen atoms (F, Cl, Br, I).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organohalogen compound, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate through atoms to find halogens bonded to carbon
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in [9, 17, 35, 53]:  # Fluorine, Chlorine, Bromine, Iodine
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:  # Carbon
                    return True, f"Halogen found bonded to carbon at atom index {atom.GetIdx()}"
    
    return False, "No carbon-halogen bonds found"