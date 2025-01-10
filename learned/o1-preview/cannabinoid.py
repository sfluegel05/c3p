"""
Classifies: CHEBI:67194 cannabinoid
"""
"""
Classifies: cannabinoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is a cannabinoid based on its SMILES string.
    Cannabinoids contain a resorcinol moiety linked to a terpenoid moiety, 
    often forming heterocyclic rings containing oxygen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cannabinoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define resorcinol moiety (benzene ring with hydroxyl groups at positions 1 and 3)
    resorcinol_smarts = "Oc1cccc(O)c1"
    resorcinol_pattern = Chem.MolFromSmarts(resorcinol_smarts)
    if not mol.HasSubstructMatch(resorcinol_pattern):
        return False, "No resorcinol moiety found"
    
    # Look for alkyl side chain attached to the benzene ring
    alkyl_chain_smarts = "c1ccc([CX4][CX4][CX4])cc1"
    alkyl_chain_pattern = Chem.MolFromSmarts(alkyl_chain_smarts)
    if not mol.HasSubstructMatch(alkyl_chain_pattern):
        return False, "No alkyl side chain attached to aromatic ring"
    
    # Look for heterocyclic ring containing oxygen (e.g., pyran ring)
    heterocycle_oxygen_smarts = "c1cc2OC(CCc2c1)"
    heterocycle_oxygen_pattern = Chem.MolFromSmarts(heterocycle_oxygen_smarts)
    if not mol.HasSubstructMatch(heterocycle_oxygen_pattern):
        return False, "No heterocyclic ring containing oxygen found"
    
    # Additional check for terpenoid side chain (alkyl chain attached to heterocycle)
    terpenoid_chain_smarts = "C[C@H](C)CC=C"
    terpenoid_chain_pattern = Chem.MolFromSmarts(terpenoid_chain_smarts)
    if not mol.HasSubstructMatch(terpenoid_chain_pattern):
        return False, "No terpenoid side chain found"

    return True, "Molecule contains structural features characteristic of cannabinoids"