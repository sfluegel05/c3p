"""
Classifies: CHEBI:65111 3-substituted propionyl-CoA(4-)
"""
"""
Classifies: 3-substituted propionyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_substituted_propionyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-substituted propionyl-CoA(4-) based on its SMILES string.
    This involves finding the Coenzyme A structure and ensuring substitution on the third carbon
    of the propionyl group.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-substituted propionyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the general CoA backbone SMARTS pattern
    coa_pattern = Chem.MolFromSmarts("C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA backbone found"

    # Define the 3-substituted propionyl pattern (R-C(=O)-S-C-C)
    # Assuming the substitution occurs on the first carbon of the chain
    propionyl_substituted_pattern = Chem.MolFromSmarts("C(=O)SCC[R]")
    if not mol.HasSubstructMatch(propionyl_substituted_pattern):
        return False, "No 3-substituted propionyl group found"

    # If both patterns match, assume it is a 3-substituted propionyl-CoA(4-)
    return True, "Contains CoA backbone with a 3-substituted propionyl group"