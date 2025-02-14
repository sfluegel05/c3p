"""
Classifies: CHEBI:65323 monoterpenoid indole alkaloid
"""
"""
Classifies: CHEBI:25680 monoterpenoid indole alkaloid
A terpenoid indole alkaloid which is biosynthesised from L-tryptophan and diisoprenoid 
(usually secolaganin) building blocks.
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_monoterpenoid_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid indole alkaloid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoterpenoid indole alkaloid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for indole substructure
    indole_pattern = Chem.MolFromSmarts("c1ccc2c(c1)nc3ccccc23")
    if not mol.HasSubstructMatch(indole_pattern):
        return False, "No indole substructure found"
    
    # Look for terpenoid substructure
    terpenoid_pattern = Chem.MolFromSmarts("C=CC(C)CC")
    if not mol.HasSubstructMatch(terpenoid_pattern):
        return False, "No terpenoid substructure found"
    
    # Check for connectivity between indole and terpenoid
    connected = any(atom.GetAtomicNum() == 7 and atom.GetDegree() > 3 for atom in mol.GetAtoms())
    if not connected:
        return False, "Indole and terpenoid not connected"
    
    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 600:
        return False, "Molecular weight outside expected range (200-600 Da)"
    
    # Count nitrogen atoms (should have 2-3)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 2 or n_count > 3:
        return False, f"Unexpected number of nitrogen atoms ({n_count})"
    
    return True, "Contains indole and terpenoid substructures connected"