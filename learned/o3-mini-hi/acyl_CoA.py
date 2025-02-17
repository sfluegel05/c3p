"""
Classifies: CHEBI:17984 acyl-CoA
"""
"""
Classifies: Acyl-CoA
Definition: A thioester that results from the formal condensation of the thiol group 
of coenzyme A with the carboxy group of any carboxylic acid.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an acyl-CoA based on its SMILES string.
    An acyl-CoA contains a thioester group (C(=O)S) and the characteristic coenzyme A fragment.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is an acyl-CoA, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for a thioester: a carbonyl (C=O) directly bound to a sulfur.
    thioester_pattern = Chem.MolFromSmarts("[CX3](=O)[SX2]")
    if thioester_pattern is None:
        return False, "Error creating thioester SMARTS pattern"
    
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester group (C(=O)S) found"

    # Define a simplified SMARTS fragment for the Coenzyme A moiety.
    # This fragment (SCCNC(=O)CCNC(=O)) is chosen to capture the typical motif found in acyl-CoA.
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    if coa_pattern is None:
        return False, "Error creating Coenzyme A SMARTS pattern"
    
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A moiety not detected based on the expected fragment"

    # Optionally, one could check if the thioester linkage is adjacent to the CoA section.
    # For this example we assume that if a thioester and the characteristic CoA fragment are present,
    # the compound falls under acyl-CoA.
    
    return True, "Contains both a thioester (C(=O)S) and a Coenzyme A moiety"

# Example usage:
# status, reason = is_acyl_CoA("[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC([C@H](C(NCCC(NCCSC(=O)/C=C/CCCCCCCCCCCCC)=O)=O)O)(C)C)(=O)O)...")
# print(status, reason)