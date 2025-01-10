"""
Classifies: CHEBI:50998 trans-2-enoyl-CoA
"""
"""
Classifies: trans-2-enoyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_trans_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a trans-2-enoyl-CoA based on its SMILES string.
    These are CoA thioesters with a trans double bond at position 2 of the acyl chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trans-2-enoyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for CoA moiety
    # Look for adenine nucleobase
    adenine_pattern = Chem.MolFromSmarts("c1nc(N)c2ncnc2n1")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "No CoA moiety found (missing adenine)"

    # Look for phosphate groups characteristic of CoA
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)(O)O")
    if len(mol.GetSubstructMatches(phosphate_pattern)) < 2:
        return False, "No CoA moiety found (insufficient phosphate groups)"

    # Look for thioester linkage
    thioester_pattern = Chem.MolFromSmarts("C(=O)SC")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"

    # Look for trans double bond adjacent to carbonyl
    # The pattern \C=C\C(=O)S captures trans configuration
    trans_enoyl_pattern = Chem.MolFromSmarts("\\C=C\\C(=O)S")
    cis_enoyl_pattern = Chem.MolFromSmarts("/C=C\\C(=O)S")
    
    has_trans = mol.HasSubstructMatch(trans_enoyl_pattern)
    has_cis = mol.HasSubstructMatch(cis_enoyl_pattern)
    
    if not (has_trans or has_cis):
        return False, "No α,β-unsaturated thioester found"
    
    if has_cis and not has_trans:
        return False, "Found cis-enoyl-CoA instead of trans"

    # Additional check for pantetheine arm of CoA
    pantetheine_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    if not mol.HasSubstructMatch(pantetheine_pattern):
        return False, "Missing pantetheine arm of CoA"

    return True, "Contains CoA moiety with trans-2-enoyl group"