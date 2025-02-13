"""
Classifies: CHEBI:23437 cyanohydrin
"""
"""
Classifies: CHEBI:35667 cyanohydrin
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_cyanohydrin(smiles: str):
    """
    Determines if a molecule is a cyanohydrin based on its SMILES string.
    A cyanohydrin is an alpha-hydroxynitrile resulting from the formal addition 
    of hydrogen cyanide to the C=O bond of an aldehyde or ketone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyanohydrin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carbon with both OH and CN groups attached
    # [C;X4] = sp3 carbon
    # [OX2H1] = hydroxy group (-OH)
    # [CX2]#[NX1] = cyano group (-C≡N)
    cyanohydrin_pattern = Chem.MolFromSmarts('[C;X4]([OX2H1])([*,#1])([*,#1])[CX2]#[NX1]')
    
    if not mol.HasSubstructMatch(cyanohydrin_pattern):
        return False, "No carbon found with both hydroxy and cyano groups attached"
    
    # Count number of such patterns - should be exactly one
    matches = mol.GetSubstructMatches(cyanohydrin_pattern)
    if len(matches) > 1:
        return False, "Multiple carbons with OH and CN groups found"
    
    # Count total number of CN groups
    cyano_pattern = Chem.MolFromSmarts('[CX2]#[NX1]')
    cyano_matches = mol.GetSubstructMatches(cyano_pattern)
    if len(cyano_matches) > 1:
        return False, "Multiple cyano groups found"
        
    # Count total number of OH groups
    hydroxy_pattern = Chem.MolFromSmarts('[OX2H1]')
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    if len(hydroxy_matches) > 1:
        return False, "Multiple hydroxy groups found"
    
    # Verify carbon hybridization
    match = matches[0]
    c_atom = mol.GetAtomWithIdx(match[0])
    if c_atom.GetHybridization() != Chem.HybridizationType.SP3:
        return False, "Carbon with OH and CN must be sp3 hybridized"
        
    # Check that molecule is not too small (should be larger than HCN)
    if mol.GetNumAtoms() < 4:
        return False, "Molecule too small to be a cyanohydrin"
        
    return True, "Contains carbon with both hydroxy and cyano groups in correct configuration"