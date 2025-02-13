"""
Classifies: CHEBI:67194 cannabinoid
"""
from rdkit import Chem

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is a cannabinoid based on its SMILES string.
    Cannabinoids are identified by complex aromatic ring systems, long carbon chains, 
    and functional groups like hydroxyls, esters, amides, or ethers.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cannabinoid, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Expanded patterns for known cannabinoid core structures
    cannabinoid_core_patterns = [
        Chem.MolFromSmarts('[cR2]1[cR2][cR2][cR2][cR2][cR2]1[C]=C'), # Aromatic ring likely with double bond
        Chem.MolFromSmarts('[C@H]1CCC(C)=C[C@H]1C'),                  # Cyclohexene with sterochemistry
        Chem.MolFromSmarts('Cc1cc(O)c(O)cc1'),                       # Phenol with additional OH or methyl groups
    ]

    # Check for cannabinoid core structures
    for pattern in cannabinoid_core_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Cannabinoid-like core structure detected"
    
    # Additional characteristic functional groups in cannabinoids
    functional_groups = [
        Chem.MolFromSmarts('[OH]'),             # Hydroxyl group
        Chem.MolFromSmarts('[OX2H0]'),          # Ether linkage
        Chem.MolFromSmarts('OC(=O)'),           # Ester (carboxylate)
        Chem.MolFromSmarts('NC(=O)'),           # Amide
        Chem.MolFromSmarts('[C]CC=C'),          # Long alkenyl or alkyl chains
    ]
    
    # Check for at least one common functional group
    for pattern in functional_groups:
        if mol.HasSubstructMatch(pattern):
            return True, "Characteristic functional group indicative of cannabinoids detected"
    
    return False, "No characteristic patterns of cannabinoids detected"