"""
Classifies: CHEBI:60245 inositol phosphoceramide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_inositol_phosphoceramide(smiles: str):
    """
    Determines if a molecule is an inositol phosphoceramide based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an inositol phosphoceramide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for charged species
    if "-" in smiles or "+" in smiles:
        return False, "Charged species not included in definition"

    # Check for inositol ring (cyclohexane with OH groups, allowing substitution)
    # More flexible pattern that allows for substitution at O positions
    inositol_pattern = Chem.MolFromSmarts("[OX2,OH1][C@H]1[C@H]([OX2,OH1])[C@H]([OX2,OH1])[C@H]([OX2,OH1])[C@H]([OX2,OH1])[C@H]1[OX2,OH1]")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring found"

    # Check for single phosphodiester bridge (-O-P(=O)(O)-O-)
    phosphodiester_pattern = Chem.MolFromSmarts("[O]-[P](=O)([O])[O]")
    phosphate_matches = len(mol.GetSubstructMatches(phosphodiester_pattern))
    if phosphate_matches == 0:
        return False, "No phosphodiester bridge found"
    if phosphate_matches > 1:
        return False, "Multiple phosphate groups found"

    # Check for ceramide core structure
    # Look for sphingoid base with specific stereochemistry
    ceramide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[CX4][CX4]")
    if not mol.HasSubstructMatch(ceramide_pattern):
        return False, "No ceramide core structure found"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Should have long carbon chains (ceramide part)
    if c_count < 20:
        return False, "Carbon count too low for ceramide chains"

    # Check for minimum oxygen count
    min_oxygens = 8  # 6 from inositol, 1 amide, 1 phosphate
    if o_count < min_oxygens:
        return False, "Not enough oxygen atoms for required functional groups"

    # Verify long alkyl chains characteristic of ceramides
    alkyl_chain_pattern = Chem.MolFromSmarts("CCCCCCCCC")  # At least 9 carbons in chain
    alkyl_chains = len(mol.GetSubstructMatches(alkyl_chain_pattern))
    if alkyl_chains < 2:
        return False, "Missing long alkyl chains characteristic of ceramides"

    # Check for mannose residue (optional)
    mannose_pattern = Chem.MolFromSmarts("O[C@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]1O")
    has_mannose = mol.HasSubstructMatch(mannose_pattern)

    # Ensure proper connectivity
    # The phosphate should connect inositol to ceramide
    connectivity_pattern = Chem.MolFromSmarts("[C@H]1[C@H][C@H][C@H][C@H][C@H]1O[P]([O])(=O)OCC[C@H]([NH]C(=O))")
    if not mol.HasSubstructMatch(connectivity_pattern):
        return False, "Incorrect connectivity between inositol and ceramide"

    base_description = "Contains inositol ring with correct stereochemistry, single phosphodiester bridge, and ceramide moiety"
    if has_mannose:
        return True, f"{base_description} with mannose residue"
    return True, base_description