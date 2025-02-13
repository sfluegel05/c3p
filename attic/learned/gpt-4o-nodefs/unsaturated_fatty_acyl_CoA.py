"""
Classifies: CHEBI:51006 unsaturated fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_unsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acyl-CoA based on its SMILES string.
    An unsaturated fatty acyl-CoA contains an unsaturated fatty acid chain attached to CoA via a thioester linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an unsaturated fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Improved CoA pattern covering more of the key structural elements
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA moiety found"

    # Thioester linkage pattern: C(=O)S
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"

    # Unsaturated fatty acid chain pattern: presence of C=C double bonds
    # Assuming at least one double bond is necessary
    unsaturated_pattern = Chem.MolFromSmarts("C=C")
    if not mol.HasSubstructMatch(unsaturated_pattern):
        return False, "No unsaturated bonds found"

    # Additional structural validation may be required to fully verify authenticity

    return True, "Molecule is an unsaturated fatty acyl-CoA"

# You can test the function with a SMILES string from the provided examples