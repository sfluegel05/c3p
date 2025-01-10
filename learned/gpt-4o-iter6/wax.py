"""
Classifies: CHEBI:73702 wax
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_wax(smiles: str):
    """
    Determines if a molecule is a wax based on its SMILES string.
    A wax typically contains long alkyl chains connected via an ester linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a wax, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
   
    # Check for ester group presence
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H0]")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester linkages found"
    
    # Identify long alkyl chains
    # Pattern should match flexible long alkyl chains (e.g., both saturated and unsaturated)
    # Allow some flexibility for branches but focus on length
    long_alkyl_pattern = Chem.MolFromSmarts("[R]CCCCCCCCCCCC[CX3](=O)[OX2H0]")
    match = mol.GetSubstructMatches(long_alkyl_pattern)
    if len(match) < 1:
        return False, "No sufficient long alkyl chains found"
    
    # Check for molecular weight indicative of waxes
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 450:  # Further adjusted threshold based on common examples
        return False, "Molecular weight too low for a typical wax compound"
    
    # Exclude structures with elements uncommon in waxes
    non_wax_elements = {7, 15, 16, 17, 35, 53}  # Include only typical elements like C, H, O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in non_wax_elements:
            return False, f"Contains elements atypical of waxes: {atom.GetSymbol()}"
    
    return True, "Contains characteristic ester linkage and long alkyl chains typical of waxes"