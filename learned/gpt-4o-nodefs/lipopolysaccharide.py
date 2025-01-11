"""
Classifies: CHEBI:16412 lipopolysaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lipopolysaccharide(smiles: str):
    """
    Determines if a molecule is a lipopolysaccharide based on its SMILES string.
    The identification is based on generic features common in lipopolysaccharides like
    long lipid chains and polysaccharide components. This is, however, a rudimentary
    test due to the complexity and diversity of LPS structures.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule shows rudimentary LPS features, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for long carbon chains (indicating lipid presence)
    carbon_chain_pattern = Chem.MolFromSmarts("[C;R0][C;R0][C;R0][C;R0][C;R0][C;R0]")
    long_chain_matches = mol.HasSubstructMatch(carbon_chain_pattern)
    
    if not long_chain_matches:
        return False, "No long carbon chains indicative of lipid present"
    
    # Check for saccharide units (looking for rings with multiple OH groups)
    saccharide_pattern = Chem.MolFromSmarts("C1O[C@H](CO)[C@@H](O)[C@@H](O)[C@@H]1O")
    saccharide_matches = mol.GetSubstructMatches(saccharide_pattern)
    
    if len(saccharide_matches) < 1:
        return False, "No polysaccharide-like structures detected"
    
    # Check molecular size - LPS structures are generally large
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 1000:
        return False, "Molecular weight is too low for typical LPS"

    return True, "Contains long lipid chains and polysaccharide-like structures indicating potential LPS"