"""
Classifies: CHEBI:65323 monoterpenoid indole alkaloid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid_indole_alkaloid(smiles: str):
    """
    Classifies if a given molecule is a monoterpenoid indole alkaloid (MIA)
    based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is likely an MIA, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for indole moiety
    indole_pattern = Chem.MolFromSmarts("c1ccc2[nH]ccc2c1")
    
    # Check for indole moiety
    if not mol.HasSubstructMatch(indole_pattern):
        return False, "No indole moiety detected"

    # Look for complex ring system: at least 3 rings
    ri = mol.GetRingInfo()
    if len(ri.AtomRings()) < 3:
        return False, f"Insufficient ring complexity, found {len(ri.AtomRings())} rings"

    # Check for tryptophan-derived amine or additional ester/ether linkages
    # This is a heuristic approach based on structural diversity
    ester_or_ether = Chem.MolFromSmarts("[CX4](=O)[OX2,CX3H]")
    if not mol.HasSubstructMatch(ester_or_ether):
        return False, "No ester or ether linkage typical for MIAs"

    # Check for a molecular complexity heuristic (e.g., 15+ heavy atoms)
    heavy_atom_count = rdMolDescriptors.CalcNumHeavyAtoms(mol)
    if heavy_atom_count < 15:
        return False, "Molecule not complex enough to be an MIA"
    
    return True, "Molecule has features consistent with monoterpenoid indole alkaloids"