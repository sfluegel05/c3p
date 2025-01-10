"""
Classifies: CHEBI:73702 wax
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_wax(smiles: str):
    """
    Determines if a molecule is a wax based on its SMILES string.
    A wax typically contains long chain hydrocarbons connected via an ester linkage,
    and is often malleable at room temperature.

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
    
    # Check for long alkyl chains with flexibility for unsaturation
    long_alkyl_pattern = Chem.MolFromSmarts("[R]CCCCCCCCCCCCCCCC[CX3](=O)[OX2H0]")  # extended chain length
    alkyl_chains = mol.GetSubstructMatches(long_alkyl_pattern)
    if len(alkyl_chains) < 1:
        return False, "No sufficient long alkyl chains found"
    
    # Add flexibility for unsaturated hydrocarbons
    unsaturated_alkyl_pattern = Chem.MolFromSmarts("C=C[-CH2-]{6,}[CX3](=O)[OX2H0]")  # Accounts for double bonds
    if not mol.HasSubstructMatch(unsaturated_alkyl_pattern) and len(alkyl_chains) < 2:
        return False, "Insufficient long alkyl or unsaturated chains matching the typical wax structure."

    # Adjusted molecular weight check for typical wax compounds
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for a typical wax compound"

    # Exclude structures with elements uncommon in waxes
    non_wax_elements = {7, 15, 16, 17, 35, 53}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in non_wax_elements:
            return False, f"Contains elements atypical of waxes: {atom.GetSymbol()}"

    return True, "Contains characteristic ester linkage and sufficient alkyl chains typical of waxes"