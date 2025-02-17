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
    carbon–fluorine (C–F) bond. To capture all possible C–F bonds (including bonds 
    in aromatic systems), we use the SMARTS pattern "[#6]~[F]". This pattern matches 
    any carbon atom (atomic number 6) connected by any bond (~) to a fluorine atom.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule contains at least one carbon–fluorine bond, False otherwise.
        str: Explanation of the classification.
    """
    # Convert the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a carbon–fluorine bond.
    # Using "~" (any bond) instead of "-" (single bond) improves matching in cases such as aromatic systems.
    cf_pattern = Chem.MolFromSmarts("[#6]~[F]")
    if cf_pattern is None:
        return False, "Failed to create SMARTS pattern for C–F bond"
    
    # Check if the molecule contains any substructure matching the C–F bond pattern.
    if mol.HasSubstructMatch(cf_pattern):
        return True, "Contains at least one carbon–fluorine bond"
    else:
        return False, "No carbon–fluorine bond found"

# Example testing code; remove or comment out in production.
if __name__ == "__main__":
    test_smiles = [
        "C[C@@H]1CN(C(=O)C2=C(C=CC(=C2)NC(=O)C3=CC=NC=C3)O[C@@H]1CN(C)CC4=CC=C(C=C4)C(F)(F)F",  # CF bonds present
        "CCc1ccc(Nc2c(F)c(F)cc(F)c2F)c(CC(O)=O)c1",  # Contains several C–F bonds
        "N1(C=C(C2=C1C=CC=C2)C(=O)C3=C(C=CC=C3)I)CCCCCF",  # Contains a fluoropentyl group
        "NC(=O)CF",  # 2-fluoroacetamide; clearly contains a C–F bond
        "ClC1=CC=C(C=2N=C(ON2)C3=C(NN=C3)C4=CC=CC=C4)C=C1",  # Does not contain a C–F bond
        # Note: For groups like –SF5 (pentafluorosulfanyl), even though the molecule contains many F atoms,
        # none of them are directly attached to carbon so the compound should not be classified as an organofluorine.
        "CC1=CC(C)=C(NC(=O)C2=CC(=CC(=C2)C#N)S(F)(F)(F)(F)F)C=C1N1C=CN2N=C(C=C12)C1=CC=CN=C1",
    ]
    for s in test_smiles:
        result, reason = is_organofluorine_compound(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")