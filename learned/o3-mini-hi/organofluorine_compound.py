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
    carbon–fluorine (C–F) bond. This function uses a SMARTS pattern to identify any bond 
    between a carbon atom and a fluorine atom.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule contains at least one carbon–fluorine bond, False otherwise.
        str: Explanation for the classification.
    """
    # Convert the SMILES string into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a carbon-fluorine bond: any carbon atom directly bound to a fluorine.
    cf_pattern = Chem.MolFromSmarts("[#6]-[F]")
    if cf_pattern is None:
        return False, "SMARTS pattern creation failed"
    
    # Use RDKit's substructure search to check for a C-F bond.
    if mol.HasSubstructMatch(cf_pattern):
        return True, "Contains at least one carbon-fluorine bond"
    else:
        return False, "No carbon-fluorine bond found"

# Example usage for testing; remove or comment out before production.
if __name__ == "__main__":
    test_smiles = [
        "C[C@@H]1CN(C(=O)C2=C(C=CC(=C2)NC(=O)C3=CC=NC=C3)O[C@@H]1CN(C)CC4=CC=C(C=C4)C(F)(F)F)[C@@H](C)CO",  # Contains C-F bond
        "CCc1ccc(Nc2c(F)c(F)cc(F)c2F)c(CC(O)=O)c1",  # Contains C-F bond
        "N1(C=C(C2=C1C=CC=C2)C(=O)C3=C(C=CC=C3)I)CCCCCF",  # Contains C-F bond (fluoropentyl group)
        "ClC1=CC=C(C=2N=C(ON2)C3=C(NN=C3)C4=CC=CC=C4)C=C1",  # Does not contain any C-F bond
        # Example with a pentafluorosulfanyl group (no C-F bond because F is bound to S not C)
        "CC1=CC(C)=C(NC(=O)C2=CC(=CC(=C2)C#N)S(F)(F)(F)(F)F)C=C1N1C=CN2N=C(C=C12)C1=CC=CN=C1",
    ]
    for s in test_smiles:
        result, reason = is_organofluorine_compound(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")