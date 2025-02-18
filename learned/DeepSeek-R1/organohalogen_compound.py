"""
Classifies: CHEBI:17792 organohalogen compound
"""
"""
Classifies: organohalogen compound (CHEBI: orgHalogen)
"""
from rdkit import Chem

def is_organohalogen_compound(smiles: str):
    """
    Determines if a molecule is an organohalogen compound based on the presence of at least one carbon-halogen bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organohalogen compound, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Atomic numbers for halogens: F(9), Cl(17), Br(35), I(53), At(85)
    halogens = {9, 17, 35, 53, 85}
    
    # Iterate through all atoms in the molecule
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Check only carbon atoms
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() in halogens:
                    return True, "Contains at least one carbon-halogen bond"
    
    return False, "No carbon-halogen bonds found"