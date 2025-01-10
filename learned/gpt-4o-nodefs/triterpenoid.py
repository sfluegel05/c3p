"""
Classifies: CHEBI:36615 triterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_triterpenoid(smiles: str):
    """
    Determines if a molecule is a triterpenoid based on its SMILES string.
    Triterpenoids typically consist of C30 units derived from isoprene, which can be cyclized and functionalized.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triterpenoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check molecular weight and carbon count
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    num_carbon = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if mol_weight < 400 or mol_weight > 1000:  # More flexibility for functionalizations
        return False, f"Molecular weight of {mol_weight} is not typical for triterpenoids"
    if num_carbon < 25 or num_carbon > 35:
        return False, f"Number of carbon atoms ({num_carbon}) not typical for triterpenoids"

    # Check for diverse ring systems (general triterpenoid structure involves multiple rings)
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count < 4 or ring_count > 6:  # Many triterpenoids are pentacyclic
        return False, f"Number of rings ({ring_count}) is not typical for triterpenoids"

    # Identify common motifs in triterpenoids using relaxed SMARTS pattern
    # Using a generic pentacyclic pattern
    pentacyclic_smarts = "C1CCC2CCCC3CCCC4CCCC5C1C2345"
    substruct_mol = Chem.MolFromSmarts(pentacyclic_smarts)
    if mol.HasSubstructMatch(substruct_mol):
        # Additional known functional groups (for diversity)
        functional_groups_smarts = ["[OH]", "[=O]", "[O;R]"]  # Hydroxyl, ketone, cyclic ethers
        for pattern in functional_groups_smarts:
            fg_mol = Chem.MolFromSmarts(pattern)
            if mol.HasSubstructMatch(fg_mol):
                return True, "Structure matches characteristics of triterpenoids (pentacyclic core + functional groups)"
        return True, "Structure matches pentacyclic core characteristics of triterpenoids but lacks notable functional groups"

    return False, "Structure doesn't match common triterpenoid motifs"