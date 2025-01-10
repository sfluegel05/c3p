"""
Classifies: CHEBI:37141 organobromine compound
"""
from rdkit import Chem

def is_organobromine_compound(smiles: str):
    """
    Determines if a molecule is an organobromine compound based on its SMILES string.
    An organobromine compound contains at least one carbon-bromine bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organobromine compound, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for at least one bromine atom bonded to a carbon atom
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 35:  # Bromine's atomic number is 35
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:  # Carbon's atomic number is 6
                    return True, "Contains at least one carbon-bromine bond"
    
    return False, "No carbon-bromine bond found"

# Example usage
smiles_example = "Brc1ccc(Oc2cc(Br)ccc2Br)c(Br)c1"
is_organobromine, reason = is_organobromine_compound(smiles_example)
print(is_organobromine, reason)