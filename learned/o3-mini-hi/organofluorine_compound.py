"""
Classifies: CHEBI:37143 organofluorine compound
"""
"""
Classifies: Organofluorine Compound
Definition: An organofluorine compound is any compound containing at least one carbon-fluorine bond.
"""

from rdkit import Chem

def is_organofluorine_compound(smiles: str):
    """
    Determines if a molecule is an organofluorine compound based on its SMILES string.
    
    An organofluorine compound is defined as any compound containing at least one carbon-fluorine (C–F) bond.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule contains at least one C–F bond, False otherwise.
        str: Reason for the classification.
    """
    # Convert the SMILES string into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Iterate over each bond in the molecule.
    # Check if at least one bond connects a carbon (atomic number 6) to a fluorine (atomic number 9) atom.
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        if ((atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 9) or 
            (atom1.GetAtomicNum() == 9 and atom2.GetAtomicNum() == 6)):
            return True, "Contains at least one carbon-fluorine bond"
    
    # If no C–F bond was found, return False with a reason.
    return False, "No carbon-fluorine bond found"

# For example usage (this part can be removed or commented out in production):
if __name__ == "__main__":
    test_smiles = [
        "NC(=O)CF",  # 2-fluoroacetamide (should be True)
        "ClC1=CC=C(C=2N=C(ON2)C3=C(NN=C3)C4=CC=CC=C4)C=C1",  # A compound without a C–F bond (should be False)
    ]
    for s in test_smiles:
        result, reason = is_organofluorine_compound(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")