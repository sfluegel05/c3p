"""
Classifies: CHEBI:36092 clavulone
"""
from rdkit import Chem

def is_clavulone(smiles: str):
    """
    Determines if a molecule is a clavulone based on its SMILES string.
    A clavulone is a class of esterified prostanoids obtained from marine corals.
    
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
    
    # Check for the presence of ester groups
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, "Insufficient ester groups found"
    
    # Check for cyclic structures possibly linked by esters
    cyclic_pattern = Chem.MolFromSmarts("C1CCCCC1")
    if not mol.HasSubstructMatch(cyclic_pattern):
        return False, "No prominent cyclic structure found"
    
    # Check for conjugated diene arrangements flexible for marine prostanoids
    diene_patterns = [
        Chem.MolFromSmarts("C=CC=C"),
        Chem.MolFromSmarts("C=CC")
    ]
    has_diene = any(mol.HasSubstructMatch(dp) for dp in diene_patterns)
    if not has_diene:
        return False, "No suitable conjugated diene arrangement found"
    
    # Check for long-chain hydrocarbons typical of marine-origin compounds
    long_chain_pattern = Chem.MolFromSmarts("C=C-C-C=C")
    has_long_chain = mol.HasSubstructMatch(long_chain_pattern)
    if not has_long_chain:
        return False, "No suitable long carbon chain found"
    
    # Recognize halogens that may be indicative of specific clavulone subtypes
    halogen_pattern = Chem.MolFromSmarts("[Cl,Br,I]")
    has_halogen = mol.HasSubstructMatch(halogen_pattern)
    
    if has_halogen:
        return True, "Matches clavulone structure with halogens"
    else:
        return True, "Matches clavulone structure possibly without halogens"

# Example debug print for iodine-containing clavulone structure
# result, reason = is_clavulone("IC1=C[C@](O)(C/C=C\\CCCCC)/C(/C1=O)=C\\C=C/CCCC(OC)=O")
# print(result, reason)