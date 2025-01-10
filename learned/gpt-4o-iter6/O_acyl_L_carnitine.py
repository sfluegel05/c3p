"""
Classifies: CHEBI:75659 O-acyl-L-carnitine
"""
from rdkit import Chem

def is_O_acyl_L_carnitine(smiles: str):
    """
    Determines if a molecule is an O-acyl-L-carnitine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an O-acyl-L-carnitine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for quaternary ammonium ion pattern [N+](C)(C)C
    quaternary_ammonium_pattern = Chem.MolFromSmarts("[N+](C)(C)C")
    if not mol.HasSubstructMatch(quaternary_ammonium_pattern):
        return False, "No quaternary ammonium ion pattern found"

    # Look for ester linkage pattern OC(=O)
    ester_linkage_pattern = Chem.MolFromSmarts("OC(=O)")
    if not mol.HasSubstructMatch(ester_linkage_pattern):
        return False, "No ester linkage pattern found"

    # Check for chiral center ensuring L-configuration
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    has_l_configuration = any((symbol == "C" and chirality in {"R", "S"}) for (idx, chirality) in chiral_centers for symbol in mol.GetAtomWithIdx(idx).GetSymbol())

    if not has_l_configuration:
        return False, "No chiral centers ensuring L-configuration found"

    return True, "Contains quaternary ammonium ion, ester linkage, and L-configuration"