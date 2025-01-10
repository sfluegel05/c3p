"""
Classifies: CHEBI:67197 endocannabinoid
"""
from rdkit import Chem

def is_endocannabinoid(smiles: str):
    """
    Determines if a molecule is an endocannabinoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an endocannabinoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for long polyunsaturated chains with more flexible pattern
    unsaturated_chain_pattern = Chem.MolFromSmarts("C=C-C=C-C=C")  # Partial structure for more chains
    if not mol.HasSubstructMatch(unsaturated_chain_pattern):
        return False, "No polyunsaturated long carbon chain typical of endocannabinoids found."
    
    # Check for key functional groups with variations
    ethanolamine_pattern = Chem.MolFromSmarts("NCCO")
    amide_pattern = Chem.MolFromSmarts("NC(=O)")
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    has_key_group = any(mol.HasSubstructMatch(pattern) for pattern in [ethanolamine_pattern, amide_pattern, glycerol_pattern])

    if not has_key_group:
        return False, "No identified endocannabinoid functional group found (ethanolamine, amide, glycerol)."

    # Check linkages, focusing on typical endocannabinoid pattern
    ether_pattern = Chem.MolFromSmarts("COC")
    ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    amide_linkage_pattern = Chem.MolFromSmarts("C(=O)N")
    
    has_linkages = any(mol.HasSubstructMatch(pattern) for pattern in [ether_pattern, ester_pattern, amide_linkage_pattern])

    if not has_linkages:
        return False, "No ether, ester, or amide linkages characteristic of endocannabinoids found."
    
    # Optional: Check molecular weight (common range typical for endocannabinoids 300-400 g/mol)
    # mol_wt = Chem.rdMolDescriptors.CalcExactMolWt(mol)
    # if not (300 <= mol_wt <= 400):
    #     return False, "Molecular weight not typical of endocannabinoids."

    return True, "Matches characteristics of known endocannabinoids."