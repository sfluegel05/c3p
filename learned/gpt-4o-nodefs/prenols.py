"""
Classifies: CHEBI:26244 prenols
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_prenols(smiles: str):
    """
    Determines if a molecule is a prenol based on its SMILES string.
    Prenols generally have a terminal alcohol and consist of isoprene units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prenol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for terminal alcohol group '-OH'
    alcohol_pattern = Chem.MolFromSmarts("[CX4][OX2H]")
    if not mol.HasSubstructMatch(alcohol_pattern):
        return False, "Missing terminal alcohol group"

    # Look for isoprene units: 5-carbon chain with alternating single and double bonds
    isoprene_pattern = Chem.MolFromSmarts("C(=C-C-C=C)")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) == 0:
        return False, "No isoprene units detected"

    # Calculate molecular weight - longer prenols should have significant molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 150:
        return False, "Molecular weight too low for a typical prenol"

    # Count any conjugated systems - prenols have long conjugated carbon chains
    conjugated_system_pattern = Chem.MolFromSmarts("C=C(-*)")
    conjugated_system_matches = mol.GetSubstructMatches(conjugated_system_pattern)
    if len(conjugated_system_matches) < 3:
        return False, "Does not have enough conjugated systems characteristic of prenols"
    
    return True, "Contains terminal alcohol group and isoprene units, typical of prenols"