"""
Classifies: CHEBI:35358 sulfonamide
"""
"""
Classifies: CHEBI:33567 sulfonamide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sulfonamide(smiles: str):
    """
    Determines if a molecule is a sulfonamide based on its SMILES string.
    A sulfonamide is an amide of a sulfonic acid RS(=O)2NR'2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sulfonamide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for sulfonic amide group S(=O)(=O)N
    sulfonamide_pattern = Chem.MolFromSmarts("S(=O)(=O)N")
    if not mol.HasSubstructMatch(sulfonamide_pattern):
        return False, "No sulfonic amide group found"
    
    # Check for at least one carbon attached to N
    carbon_attached_to_N = Chem.MolFromSmarts("[N;$(NS(=O)=O)]~[C]")
    if not mol.HasSubstructMatch(carbon_attached_to_N):
        return False, "No carbon attached to sulfonamide nitrogen"
    
    # Count sulfonamide groups
    sulfonamide_count = len(mol.GetSubstructMatches(sulfonamide_pattern))
    if sulfonamide_count > 1:
        return False, f"Found {sulfonamide_count} sulfonamide groups, expected 1"
    
    # Verify molecular weight range (typically 150-500 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 150 or mol_wt > 500:
        return False, f"Molecular weight {mol_wt:.2f} outside typical sulfonamide range"
    
    return True, "Contains sulfonic amide group RS(=O)2NR'2"