"""
Classifies: CHEBI:64611 ether lipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_ether_lipid(smiles: str):
    """
    Determines if a molecule is an ether lipid based on its SMILES string.
    An ether lipid is a glycerolipid where one or more carbon atoms on glycerol
    are bonded to an alkyl chain via an ether linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ether lipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Glycerol backbone pattern (O-CH2-CH-CH2-O)
    glycerol_pattern = Chem.MolFromSmarts("[OX2]C([CH2X4])([CHX4])[OX2]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Look for ether linkages (-O-C-)
    ether_linkage_pattern = Chem.MolFromSmarts("[OX2][CX4]")
    if not mol.HasSubstructMatch(ether_linkage_pattern):
        return False, "No ether linkage found"

    # Count occurences of ester linkages as well (-C(=O)-O-)
    ester_linkage_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2]")
    ester_count = len(mol.GetSubstructMatches(ester_linkage_pattern))

    # Must have at least one ether linkage
    if mol.GetSubstructMatches(ether_linkage_pattern):
        return True, f"Contains ether linkages; Ester linkage count: {ester_count}"
    else:
        return False, "Lacks ether linkage characteristic of ether lipids"

# Example usage
example_smiles = "P(OCC[N+](C)(C)C)(OC[C@H](OC(=O)C=C)COCCCCCCCCCCCCCCCCCC)([O-])=O"
result, reason = is_ether_lipid(example_smiles)
print(result, reason)