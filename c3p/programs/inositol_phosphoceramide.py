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

    # Check for inositol ring (cyclohexane with 6 oxygens attached)
    inositol_pattern = Chem.MolFromSmarts("[OH1][C@H]1[C@H]([OH1])[C@H]([OH1])[C@H]([OH1])[C@H]([OH1])[C@H]1[OH1,O]")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol ring found"

    # Check for phosphodiester bridge (-O-P(=O)(O)-O-)
    phosphodiester_pattern = Chem.MolFromSmarts("[O]-[P](=O)([O])[O]")
    if not mol.HasSubstructMatch(phosphodiester_pattern):
        return False, "No phosphodiester bridge found"

    # Check for amide bond (ceramide linkage)
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide bond (ceramide linkage) found"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Should have long carbon chains (ceramide part)
    if c_count < 20:
        return False, "Carbon count too low for ceramide chains"

    # Check for hydroxyl groups (ceramide typically has several)
    if o_count < 8:  # minimum: 6 from inositol, 1 amide, 1 phosphate
        return False, "Not enough oxygen atoms for required functional groups"

    # Optional: Check for mannose residue
    mannose_pattern = Chem.MolFromSmarts("O[C@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]1O")
    has_mannose = mol.HasSubstructMatch(mannose_pattern)

    # Verify long alkyl chains characteristic of ceramides
    alkyl_chain_pattern = Chem.MolFromSmarts("CCCCCCCC")  # At least 8 carbons in chain
    if len(mol.GetSubstructMatches(alkyl_chain_pattern)) < 2:
        return False, "Missing long alkyl chains characteristic of ceramides"

    # Count rotatable bonds to verify flexibility of chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 15:
        return False, "Not enough rotatable bonds for ceramide chains"

    base_description = "Contains inositol ring, phosphodiester bridge, and ceramide moiety"
    if has_mannose:
        return True, f"{base_description} with mannose residue"
    else:
        return True, base_description