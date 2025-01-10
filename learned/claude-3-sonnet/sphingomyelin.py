"""
Classifies: CHEBI:64583 sphingomyelin
"""
"""
Classifies: sphingomyelin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sphingomyelin(smiles: str):
    """
    Determines if a molecule is a sphingomyelin based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a sphingomyelin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phosphocholine group
    phosphocholine = Chem.MolFromSmarts('[O,OH1]-P(=O)([O-])OCC[N+](C)(C)C')
    if not mol.HasSubstructMatch(phosphocholine):
        return False, "Missing phosphocholine group"

    # Check for amide linkage
    amide = Chem.MolFromSmarts('[NX3][CX3](=O)')
    if not mol.HasSubstructMatch(amide):
        return False, "Missing amide linkage"

    # Multiple patterns for sphingoid base to account for different variants
    sphingoid_patterns = [
        # Basic sphinganine backbone
        '[CH2X4]-[CH1X4]([NX3])-[CH1X4]([OX2])-[CH2X4]',
        # Sphingosine backbone (with double bond)
        '[CH2X4]-[CH1X4]([NX3])-[CH1X4]([OX2])-[CH1]=[CH1]',
        # More general pattern for modified sphingoid bases
        '[CH2X4]-[CX4]([NX3])-[CX4]([OX2])-[#6]'
    ]
    
    found_sphingoid = False
    for pattern in sphingoid_patterns:
        sphinx = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(sphinx):
            found_sphingoid = True
            break
    
    if not found_sphingoid:
        return False, "Missing sphingoid base structure"

    # Verify phosphocholine is connected to sphingoid base
    # Look for O-CH2-CH(NH)-CH(O)-CH2-O-P pattern
    phospho_connection = Chem.MolFromSmarts('[OX2]-[CH2X4]-[CH1X4]([NX3])-[CH1X4]([OX2])-[CH2X4]-[OX2]-P')
    if not mol.HasSubstructMatch(phospho_connection):
        return False, "Phosphocholine not properly connected to sphingoid base"

    # Count key atoms
    atom_counts = {}
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        atom_counts[symbol] = atom_counts.get(symbol, 0) + 1

    # Basic composition requirements
    if atom_counts.get('N', 0) != 2:  # One from amide, one from phosphocholine
        return False, "Must have exactly 2 nitrogen atoms"
    if atom_counts.get('P', 0) != 1:
        return False, "Must have exactly 1 phosphorus atom"
    if atom_counts.get('O', 0) < 5:  # Phosphate, amide, hydroxyl
        return False, "Insufficient oxygen atoms"
    
    # Check for long carbon chains (fatty acid + sphingoid base)
    n_carbons = atom_counts.get('C', 0)
    if n_carbons < 20:
        return False, "Insufficient carbon atoms for sphingomyelin"

    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 600:
        return False, "Molecular weight too low for sphingomyelin"

    # Check for fatty acid chain
    fatty_acid = Chem.MolFromSmarts('C(=O)-[CH2][CH2][CH2]')
    if not mol.HasSubstructMatch(fatty_acid):
        return False, "Missing fatty acid chain"

    return True, "Contains sphingoid base with amide-linked fatty acid and phosphocholine group"