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
    (Improvement: Using a SMARTS query to detect C-F bonds via any connectivity type, including aromatic bonds.)

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
    
    # Define a SMARTS pattern that looks for any carbon (atomic number 6)
    # connected (by any bond type ~) to a fluorine (atomic number 9).
    cf_smarts = "[#6]~[#9]"
    cf_pattern = Chem.MolFromSmarts(cf_smarts)
    
    if mol.HasSubstructMatch(cf_pattern):
        return True, "Contains at least one carbon-fluorine bond"
    else:
        return False, "No carbon-fluorine bond found"

# Example usage (for testing; remove or comment out for production):
if __name__ == "__main__":
    test_smiles = [
        "NC(=O)CF",  # 2-fluoroacetamide, should be True
        "ClC1=CC=C(C=2N=C(ON2)C3=C(NN=C3)C4=CC=CC=C4)C=C1",  # A compound with no C–F bond, should be False
        "CCc1ccc(Nc2c(F)c(F)cc(F)c2F)c(CC(O)=O)c1",  # robenacoxib, should be True
    ]
    for s in test_smiles:
        result, reason = is_organofluorine_compound(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")