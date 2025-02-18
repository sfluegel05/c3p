"""
Classifies: CHEBI:36092 clavulone
"""
"""
Classifies: clavulone (esterified prostanoids from marine corals)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_clavulone(smiles: str):
    """
    Determines if a molecule is a clavulone based on its SMILES string.
    Clavulones are esterified prostanoids characterized by a cyclopentenone core,
    ester groups, and unsaturated carbon chains.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a clavulone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for cyclopentenone core (cyclopent-2-en-1-one)
    cyclopentenone = Chem.MolFromSmarts('C1C(=O)C=CC1')
    core_matches = mol.GetSubstructMatches(cyclopentenone)
    if not core_matches:
        return False, "No cyclopentenone core found"
    
    # Check for at least one ester group (O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts('[OX2][CX3](=[OX1])')
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester groups found"
    
    # Check for unsaturated chain (at least 4 carbons with a double bond)
    # SMARTS pattern for chain of 4+ carbons with at least one double bond
    chain_pattern = Chem.MolFromSmarts('C=CCCC')
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    if not chain_matches:
        # Alternative pattern for different double bond positions
        chain_pattern2 = Chem.MolFromSmarts('CC=CCC')
        chain_matches2 = mol.GetSubstructMatches(chain_pattern2)
        if not chain_matches2:
            return False, "No unsaturated chain of at least 4 carbons found"
    
    # Additional check for molecular weight (clavulones are typically medium-large)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for clavulone"
    
    return True, "Contains cyclopentenone core, ester groups, and unsaturated chain"