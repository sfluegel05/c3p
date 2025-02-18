"""
Classifies: CHEBI:37142 organoiodine compound
"""
"""
Classifies: CHEBI:? organoiodine compound
An organoiodine compound is defined as a compound containing at least one carbon-iodine bond.
We refine this by requiring that the iodine atom is monovalent (degree = 1) and its only neighbor is a carbon atom.
"""

from rdkit import Chem

def is_organoiodine_compound(smiles: str):
    """
    Determines if a molecule is an organoiodine compound based on its SMILES string.
    An organoiodine compound is defined as one that contains at least one carbon-iodine bond,
    where iodine is monovalently attached (degree == 1) and directly bonded to a carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as an organoiodine compound, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Iterate over all atoms in the molecule.
    for atom in mol.GetAtoms():
        # Check if the atom is iodine (atomic number 53)
        if atom.GetAtomicNum() == 53:
            # We require that the iodine atom has exactly one neighbor (monovalent iodine)
            if atom.GetDegree() == 1:
                # Get the single neighbor.
                neighbor = atom.GetNeighbors()[0]
                # Check if the neighbor is a carbon (atomic number 6)
                if neighbor.GetAtomicNum() == 6:
                    return True, "Compound contains at least one carbon-iodine bond."
    
    # If no iodine atom fulfills the criteria, we conclude the molecule is not of the class.
    return False, "Compound does not contain any carbon-iodine bonds."

# Example usage (uncomment for testing):
# compounds = [
#     "IC1=CC(F)=C(NC(=O)CC2(CCCCC2)CC(O)=O)",  # true positive
#     "C(O)C1OC(OC1)CI"                         # should be false positive if iodine is not monovalently attached
# ]
# for s in compounds:
#     result, reason = is_organoiodine_compound(s)
#     print(s, result, reason)