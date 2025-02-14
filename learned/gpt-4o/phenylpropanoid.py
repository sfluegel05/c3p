"""
Classifies: CHEBI:26004 phenylpropanoid
"""
from rdkit import Chem

def is_phenylpropanoid(smiles: str):
    """
    Determines if a molecule belongs to the phenylpropanoid class based on its SMILES string.
    The class is determined by the presence of an aromatic system, key structural linkages, and characteristic functional groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a phenylpropanoid, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Detect aromatic rings
    aromatic_ring_pattern = Chem.MolFromSmarts("c1ccccc1")  # Aromatic phenyl group
    if not mol.HasSubstructMatch(aromatic_ring_pattern):
        return False, "No aromatic system found"

    # Complex phenylpropanoid structures often have C6-C3 backbones with variations
    phenylpropanoid_varied_topology = Chem.MolFromSmarts("c1ccc(cc1)C1CC1")  # Expanding topology recognition
    if not mol.HasSubstructMatch(phenylpropanoid_varied_topology):
        return False, "No phenylpropanoid-related structure found"

    # Functional groups common in phenylpropanoids
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")  # Hydroxyl (-OH) groups
    methoxy_ether_pattern = Chem.MolFromSmarts("COc")  # Methoxy (-OCH3) typically in ethers
    ester_pattern = Chem.MolFromSmarts("C(=O)O")  # Ester groups
    olefinic_attachment_pattern = Chem.MolFromSmarts("C=C")  # Double bonds
    
    has_functional_group = any(
        mol.HasSubstructMatch(pattern) 
        for pattern in [hydroxyl_pattern, methoxy_ether_pattern, ester_pattern, olefinic_attachment_pattern]
    )
    
    if not has_functional_group:
        return False, "Characteristic functional groups for phenylpropanoids are missing"

    return True, "Contains aromatic system and expected phenylpropanoid characteristics"