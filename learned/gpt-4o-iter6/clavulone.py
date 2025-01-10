"""
Classifies: CHEBI:36092 clavulone
"""
from rdkit import Chem

def is_clavulone(smiles: str):
    """
    Determines if a molecule is a clavulone based on its SMILES string.
    A clavulone is derived from marine corals and is classified as an esterified prostanoid.
    
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
    
    # Look for any ester group pattern
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester group found"
    
    # Check for presence of cyclic esterified structures typical for prostanoids
    cyclic_esters_pattern = Chem.MolFromSmarts("C1=CC=CCC=C1C(=O)OC")
    if not mol.HasSubstructMatch(cyclic_esters_pattern):
        return False, "No cyclic ester pattern found"
    
    # Look for more flexible conjugated diene arrangements
    diene_patterns = [
        Chem.MolFromSmarts("C=C-C=C"),
        Chem.MolFromSmarts("C=C-C-C=C"),
        Chem.MolFromSmarts("C=C-C(=C)")
    ]
    has_diene = any(mol.HasSubstructMatch(dp) for dp in diene_patterns)
    if not has_diene:
        return False, "No suitable conjugated diene arrangement found"
    
    # Check long carbon chains with branching (common in marine-derived entities)
    chain_patterns = [
        Chem.MolFromSmarts("CCCCCC=CCCC"),  # Basic pattern
        Chem.MolFromSmarts("C-C-C-C=C-C-C"),  # Adjusted for branches
    ]
    has_long_chain = any(mol.HasSubstructMatch(cp) for cp in chain_patterns)
    if not has_long_chain:
        return False, "No appropriately patterned long carbon chain found"
    
    # Check for halogens indicating clavulones with halogen variance
    halogen_pattern = Chem.MolFromSmarts("[Cl,Br,I]")
    has_halogen = mol.HasSubstructMatch(halogen_pattern)
    
    if has_halogen:
        return True, "Matches clavulone structure with halogen"
    else:
        return True, "Matches clavulone structure without halogen"

# Example debug print for iodine-containing clavulone structure
# result, reason = is_clavulone("IC1=C[C@](O)(C/C=C\\CCCCC)/C(/C1=O)=C\\C=C/CCCC(OC)=O")
# print(result, reason)