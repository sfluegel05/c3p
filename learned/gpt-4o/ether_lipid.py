"""
Classifies: CHEBI:64611 ether lipid
"""
from rdkit import Chem

def is_ether_lipid(smiles: str):
    """
    Determines if a molecule is an ether lipid based on its SMILES string.
    An ether lipid features one or more ether linkages, typically replacing ester linkages in a glycerol backbone.

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

    # Detailed pattern for ether-linked glycerol backbone
    ether_glycerol_pattern = Chem.MolFromSmarts("O[C@H](COC)CO")  # Captures ether bond on glycerol
    headgroup_pattern = Chem.MolFromSmarts("C[N+](C)(C)C")  # Phosphocholine example, common in ether lipids

    # Check for ether linkage in glycerol backbone
    has_ether_glycerol = mol.HasSubstructMatch(ether_glycerol_pattern)
    
    # Check for common phospholipid headgroups
    has_headgroup = mol.HasSubstructMatch(headgroup_pattern)

    if has_ether_glycerol and has_headgroup:
        return True, "Contains ether linkage with glycerol backbone and phospholipid headgroup"
    
    return False, "No ether-linked glycerol structure or missing characteristic headgroup"

# Example usage
example_smiles = "P(OCC[N+](C)(C)C)(OC[C@H](OC(=O)C=C)COCCCCCCCCCCCCCCCCCC)([O-])=O"
result, reason = is_ether_lipid(example_smiles)
print(result, reason)