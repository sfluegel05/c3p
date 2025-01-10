"""
Classifies: CHEBI:17984 acyl-CoA
"""
from rdkit import Chem

def is_acyl_CoA(smiles: str) -> (bool, str):
    """
    Classifies the molecule as acyl-CoA based on SMILES string.
    Acyl-CoA is characterized by a fatty acid chain bound to Coenzyme A via a thioester bond.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as acyl-CoA, False otherwise.
        str: Reason for the classification.
    """

    # Parse the SMILES string to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern to identify the fundamental Coenzyme A structure
    # CoA comprises a phosphopantetheine group and ADP
    phosphopantetheine_pattern = Chem.MolFromSmarts("P(=O)(O)OP(=O)(O)OC[C@H]1O[C@H](CO)COP(=O)(O)OP(=O)(O)OC[C@H]_4O[C@H]([C@H](O)[C@@H]_4OP(=O)(O)O)n1c2nc[nH]c2n1")
    
    # Adjusted to match fully expanded CoA moieties. Patterns leftover from simplification:
    adp_pattern = Chem.MolFromSmarts("N1C=NC2=C1C(=O)NC=N2")
    coa_core_pattern = Chem.MolFromSmarts("NC(=O)CCSC(=O)")
    
    # Combine checks for the updated components
    if not mol.HasSubstructMatch(phosphopantetheine_pattern):
        return False, "CoA moiety (phosphopantetheine group) not identified"
    
    if not mol.HasSubstructMatch(adp_pattern):
        return False, "CoA moiety (adenine group) not identified"
    
    if not mol.HasSubstructMatch(coa_core_pattern):
        return False, "Core thioester-linked CoA structure not identified"

    # To confirm, searching for thioester linkage
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCCN")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester bond with CoA moiety"

    return True, "Contains CoA moiety with a thioester bond, characteristic of acyl-CoA"