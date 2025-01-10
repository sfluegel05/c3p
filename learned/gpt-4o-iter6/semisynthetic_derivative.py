"""
Classifies: CHEBI:72588 semisynthetic derivative
"""
from rdkit import Chem

def is_semisynthetic_derivative(smiles: str):
    """
    Determines if a molecule is a semisynthetic derivative based on its SMILES string.
    Heuristic detection of synthetic modifications, focusing on unique structural alerts.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule may be a semisynthetic derivative, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Expanded patterns that may indicate synthetic modifications
    # Look for amides linked to unusual groups, exotic rings, or complex modifications
    amide_pattern = Chem.MolFromSmarts("C(=O)N")  # Simple amide group
    exotic_ring_pattern = Chem.MolFromSmarts("c1[cR2][cR2][cR2][cR2][cR2]c1")  # Aromatic ring with potential modifications
    ether_pattern = Chem.MolFromSmarts("C-O-C(O)-C")  # Ethers attached to unconventional oxygens/groups

    # Check for synthetic signatures
    has_amide = mol.HasSubstructMatch(amide_pattern)
    has_exotic_ring = mol.HasSubstructMatch(exotic_ring_pattern)
    has_specific_ether = mol.HasSubstructMatch(ether_pattern)

    # Combine indicators for a more robust system
    if has_amide and (has_specific_ether or has_exotic_ring):
        return True, "Synthetic markers detected including complex amide or exotic modifications"

    # Without explicit synthesis pathway info, classification remains speculative
    return False, "No sufficient synthetic indicators in SMILES"

# Examples to test the function
print(is_semisynthetic_derivative("CO[C@H]1C[C@H](O[C@H]2[C@H](C)[C@@H](O[C@@H]3O[C@H](C)C[C@@H]([C@H]3OC(C)=O)N(C)C)[C@@H](C)C[C@@]3(CO3)C(=O)[C@H](C)[C@@H](OC(C)=O)[C@@H](C)[C@@H](C)OC(=O)[C@@H]2C)O[C@@H](C)[C@@H]1OC(C)=O"))  # Troleandomycin
print(is_semisynthetic_derivative("CCCCn1c2cc(OC)ccc2c2ccnc(C)c12"))  # N-butylharmine