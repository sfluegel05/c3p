"""
Classifies: CHEBI:17761 ceramide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ceramide(smiles: str):
    """
    Determines if a molecule is a ceramide based on its SMILES string.
    
    A ceramide is a sphingoid base derivative with an amide-linked fatty acid.
    The fatty acids are typically saturated or monounsaturated with chain lengths
    from 14 to 26 carbon atoms, and the presence of a hydroxyl group on carbon 2
    is fairly common. Ceramides are generally precursors of more complex sphingolipids.

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
    
    # Look for sphingosine backbone pattern
    sphingosine_pattern = Chem.MolFromSmarts("[CH2X4][NX3][CH2X4][CH2X4][CH2X4][CH0X3]([OH1])[CH2X4]")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No sphingosine backbone found"
    
    # Check for fatty acid chain
    fatty_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[NX3][CX4]([CH1X4])[CH2X4][CH2X4][CH2X4][CH2X4]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if not fatty_acid_matches:
        return False, "No fatty acid chain found"
    
    # Check fatty acid chain length (14-26 carbons)
    fatty_acid_chain_length = sum(1 for atom in mol.GetAtomWithIdx(fatty_acid_matches[0][2]).GetNeighbors() if atom.GetAtomicNum() == 6)
    if not (14 <= fatty_acid_chain_length <= 26):
        return False, f"Fatty acid chain length ({fatty_acid_chain_length}) outside typical range (14-26)"
    
    # Check for hydroxyl group on C2 of fatty acid chain
    c2_fatty_acid = mol.GetAtomWithIdx(fatty_acid_matches[0][3])
    has_c2_hydroxyl = any(atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == 0 for atom in c2_fatty_acid.GetNeighbors())
    
    # Check molecular weight (typically > 500 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for ceramide"
    
    # Additional checks or filtering criteria can be added here
    
    return True, "Contains sphingosine backbone with amide-linked fatty acid chain"