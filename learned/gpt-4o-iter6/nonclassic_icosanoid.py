"""
Classifies: CHEBI:61703 nonclassic icosanoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nonclassic_icosanoid(smiles: str):
    """
    Determines if a molecule is a nonclassic icosanoid based on its SMILES string.
    These are C20 fatty acid derivatives with oxygenation, excluding leukotrienes and prostanoids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nonclassic icosanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for a continuous carbon chain with approximately 20 carbons
    c_chain_pattern = Chem.MolFromSmarts("C" * 20)
    if not mol.HasSubstructMatch(c_chain_pattern):
        return False, "Does not have a C20 carbon backbone"

    # Check for oxygenation: hydroxyl, epoxy, or carboxyl groups
    oxygen_groups = ["[OH]", "[C@H]1O[C@H]1", "[O;X2]=C"]
    oxygen_count = sum(mol.HasSubstructMatch(Chem.MolFromSmarts(pat)) for pat in oxygen_groups)
    if oxygen_count < 1:
        return False, f"Insufficient oxygenation features, found {oxygen_count} oxygen groups"

    # Calculate number of carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (18 <= c_count <= 22):  # Allow a little flexibility around C20
        return False, f"Carbon count out of range: {c_count}"

    # Check for carboxyl group indicating fatty acid
    carboxylic_pattern = Chem.MolFromSmarts("C(=O)O")
    has_carboxyl_group = mol.HasSubstructMatch(carboxylic_pattern)
    if not has_carboxyl_group:
        return False, "Missing carboxyl group, typical of fatty acids"

    # Exclude typical leukotriene pattern if identifiable
    leukotriene_pattern = Chem.MolFromSmarts("CCCCC(C=CC=CC=CC=CC=CC=CC=C)C")
    if mol.HasSubstructMatch(leukotriene_pattern):
        return False, "Structure matches that of typical leukotrienes"

    # Exclude prostanoid-like structure
    prostanoid_pattern = Chem.MolFromSmarts("CC(C)CCC1C2CCC3C1C(CCC3C2)C(O)=O")
    if mol.HasSubstructMatch(prostanoid_pattern):
        return False, "Structure matches that of typical prostanoids"

    return True, "Contains characteristics of nonclassic icosanoids: C20 backbone with oxygenation"