"""
Classifies: CHEBI:17792 organohalogen compound
"""
"""
Classifies: Organohalogen Compound
A compound containing at least one carbon-halogen bond (where X is a halogen atom).
"""

from rdkit import Chem

def is_organohalogen_compound(smiles: str):
    """
    Determines if a molecule is an organohalogen compound based on its SMILES string.
    An organohalogen compound must contain at least one carbon-halogen (C-X) bond, 
    where X is one of the halogen atoms: F, Cl, Br, or I.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organohalogen compound, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the atomic numbers for carbon and chosen halogens
    carbon_atomic_num = 6
    halogen_atomic_nums = {9, 17, 35, 53}  # F, Cl, Br, I respectively
    
    # Iterate through all bonds to check for a carbon-halogen bond
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        # Check if one atom is carbon and the other is a halogen
        if (atom1.GetAtomicNum() == carbon_atomic_num and atom2.GetAtomicNum() in halogen_atomic_nums) or \
           (atom2.GetAtomicNum() == carbon_atomic_num and atom1.GetAtomicNum() in halogen_atomic_nums):
            return True, "Molecule contains a carbon-halogen bond"
    
    return False, "No carbon-halogen bond found"