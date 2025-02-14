"""
Classifies: CHEBI:64611 ether lipid
"""
from rdkit import Chem

def is_ether_lipid(smiles: str):
    """
    Determines if a molecule is an ether lipid based on its SMILES string.
    An ether lipid features one or more ether linkages replacing the usual ester linkages with a glycerol backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ether lipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string to RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the ether linkage pattern (R-O-R') that should replace ester linkage (R-C(=O)-O-R)
    ether_linkage_pattern = Chem.MolFromSmarts("COC")  # Defines a simple ether linkage from glycerol to alkyl chain

    # Define the glycerol ether pattern to ensure it's part of the backbone
    # This should capture more varied structures with different head groups after the glycerol moiety
    glycerol_ether_pattern = Chem.MolFromSmarts("C(CO)CO")  # Simplified glycerol backbone, allowing ether linkage

    # Check for presence of ether linkage
    has_ether_linkage = mol.HasSubstructMatch(ether_linkage_pattern)
    
    # Check if the molecule has a glycerol backbone element with potential ether links
    has_glycerol_structure = mol.HasSubstructMatch(glycerol_ether_pattern)
    
    # Determine if the molecule is an ether lipid
    if has_ether_linkage and has_glycerol_structure:
        return True, "Contains ether linkage with glycerol backbone"
    
    return False, "No ether-linked glycerol backbone found"

# Example usage
example_smiles = "P(OCC[N+](C)(C)C)(OC[C@H](OC(=O)C=C)COCCCCCCCCCCCCCCCCCC)([O-])=O"
result, reason = is_ether_lipid(example_smiles)
print(result, reason)