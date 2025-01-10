"""
Classifies: CHEBI:17761 ceramide
"""
"""
Classifies: CHEBI:17761 ceramide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_ceramide(smiles: str):
    """
    Determines if a molecule is a ceramide based on its SMILES string.
    A ceramide is a sphingoid base with an amide-linked fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ceramide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # More flexible sphingoid base pattern
    # Matches long carbon chain with at least one hydroxyl and an amine group
    sphingoid_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CHX4]([OH])[NH]")
    if not mol.HasSubstructMatch(sphingoid_pattern):
        # Try alternative pattern with double bonds
        sphingoid_pattern = Chem.MolFromSmarts("[CH2X4][CHX4]/[CHX4]=[CHX4]/[CHX4]([OH])[NH]")
        if not mol.HasSubstructMatch(sphingoid_pattern):
            # Try another alternative pattern with additional hydroxyl groups
            sphingoid_pattern = Chem.MolFromSmarts("[CH2X4][CHX4]([OH])[CHX4]([OH])[NH]")
            if not mol.HasSubstructMatch(sphingoid_pattern):
                return False, "No sphingoid base found"

    # Look for amide group (-C(=O)-N-)
    amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3]")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) != 1:
        return False, f"Found {len(amide_matches)} amide groups, need exactly 1"

    # Check for fatty acid chain (long carbon chain attached to amide)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") 
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 1:
        return False, f"Missing fatty acid chain, got {len(fatty_acid_matches)}"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Chains too short to be fatty acids"

    # Check molecular weight - ceramides typically >400 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for ceramide"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 20:
        return False, "Too few carbons for ceramide"
    if o_count < 2:
        return False, "Must have at least 2 oxygens (amide and hydroxyl groups)"

    return True, "Contains sphingoid base with amide-linked fatty acid chain"