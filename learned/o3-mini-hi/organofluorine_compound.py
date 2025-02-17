"""
Classifies: CHEBI:37143 organofluorine compound
"""
"""
Classifies: Organofluorine Compound
Definition: An organofluorine compound is any compound containing at least one carbon–fluorine bond.
"""

from rdkit import Chem

def is_organofluorine_compound(smiles: str):
    """
    Determines if a molecule is an organofluorine compound based on its SMILES string.

    An organofluorine compound is defined as any compound containing at least one 
    carbon–fluorine (C–F) bond. This function parses the SMILES and then searches
    through all fluorine atoms to check whether any one is directly bonded to a carbon atom.
    
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
    
    # Iterate over all atoms in the molecule.
    for atom in mol.GetAtoms():
        # Check if the atom is a fluorine (atomic number 9).
        if atom.GetAtomicNum() == 9:  
            # Loop over neighboring atoms of the fluorine.
            for neighbor in atom.GetNeighbors():
                # Check if neighbor is a carbon (atomic number 6).
                if neighbor.GetAtomicNum() == 6:
                    return True, "Contains at least one carbon-fluorine bond"
    
    # If no fluorine atom is bonded to a carbon, then it's not an organofluorine compound.
    return False, "No carbon-fluorine bond found"

# Example usage (for testing; remove or comment out for production):
if __name__ == "__main__":
    test_smiles = [
        "NC(=O)CF",  # 2-fluoroacetamide, should be True
        "CCc1ccc(Nc2c(F)c(F)cc(F)c2F)c(CC(O)=O)c1",  # robenacoxib, should be True
        "N1(C=C(C2=C1C=CC=C2)C(=O)C3=C(C=CC=C3)I)CCCCCF",  # 1-(5-fluoropentyl)-3-(2-iodobenzoyl)indole, should be True
        "ClC1=CC=C(C=2N=C(ON2)C3=C(NN=C3)C4=CC=CC=C4)C=C1",  # A compound with no C–F bond, should be False
        # Compound with pentafluorosulfanyl group (S(F)(F)(F)(F)F) but no C–F bond:
        "CC1=CC(C)=C(NC(=O)C2=CC(=CC(=C2)C#N)S(F)(F)(F)(F)F)C=C1N1C=CN2N=C(C=C12)C1=CC=CN=C1",
    ]
    for s in test_smiles:
        result, reason = is_organofluorine_compound(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")