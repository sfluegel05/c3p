"""
Classifies: CHEBI:37141 organobromine compound
"""
"""
Classifies: Organobromine compound
Definition: A compound containing at least one carbon–bromine bond.
This version inspects each bond in the molecule directly to see whether one atom is a carbon (atomic number 6)
and the other is a bromine (atomic number 35). This avoids potential pitfalls using SMARTS patterns.
"""

from rdkit import Chem

def is_organobromine_compound(smiles: str):
    """
    Determines if a molecule is an organobromine compound based on its SMILES string.
    An organobromine compound is defined as a compound containing at least one carbon–bromine bond.
    
    This function parses the molecule and then iterates over all bonds.
    If any bond is found connecting a carbon atom (atomic number 6) to a bromine atom (atomic number 35),
    the molecule is classified as an organobromine compound.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is an organobromine compound, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Iterate over every bond in the molecule
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        a1_atomic = atom1.GetAtomicNum()
        a2_atomic = atom2.GetAtomicNum()
        # Check if one atom is carbon (6) and the other is bromine (35).
        if (a1_atomic == 6 and a2_atomic == 35) or (a1_atomic == 35 and a2_atomic == 6):
            return True, "Contains at least one carbon–bromine bond"
    
    # No carbon–bromine bond found.
    return False, "No carbon–bromine bond found"