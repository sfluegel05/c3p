"""
Classifies: CHEBI:32957 lysophosphatidic acids
"""
"""
Classifies: CHEBI:16480 lysophosphatidic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_lysophosphatidic_acids(smiles: str):
    """
    Determines if a molecule is a lysophosphatidic acid based on its SMILES string.
    A lysophosphatidic acid has a glycerol backbone with one fatty acid chain attached via an ester bond and a phosphate group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lysophosphatidic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with hydroxyl groups)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
        
    # Look for exactly one ester group (-O-C(=O)-) attached to the glycerol backbone
    ester_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4](-[OX2][CX3](=[OX1]))")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups attached to glycerol backbone, need exactly 1"

    # Look for a phosphate group (P(=O)(O)(O)- or P(O)(O)=O) attached to the glycerol backbone
    phosphate_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4](-[PX4](=[OX1])([OX2])[OX2])")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) != 1:
        return False, f"Found {len(phosphate_matches)} phosphate groups attached to glycerol backbone, need exactly 1"

    # Check for fatty acid chain (long carbon chain attached to ester)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") 
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 1:
        return False, "Missing fatty acid chain"

    # Count rotatable bonds to verify long chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Chain too short to be a fatty acid"

    # Check molecular weight - LPAs typically >400 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for lysophosphatidic acid"

    return True, "Contains glycerol backbone with one fatty acid chain and a phosphate group"