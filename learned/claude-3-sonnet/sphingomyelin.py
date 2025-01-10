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

    # Check for phosphocholine group - simplified pattern
    phosphocholine = Chem.MolFromSmarts('OP(=O)([O-])OCC[N+](C)(C)C')
    if not mol.HasSubstructMatch(phosphocholine):
        return False, "Missing phosphocholine group"

    # Check for amide linkage
    amide = Chem.MolFromSmarts('[NX3][CX3](=O)')
    if not mol.HasSubstructMatch(amide):
        return False, "Missing amide linkage"

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
    if atom_counts.get('O', 0) < 5:  # At least 5 oxygens (phosphate, amide, hydroxyl)
        return False, "Insufficient oxygen atoms"
    if atom_counts.get('C', 0) < 20:  # Sphingomyelins typically have many carbons
        return False, "Insufficient carbon atoms"

    # Check for hydroxyl group(s)
    hydroxyl = Chem.MolFromSmarts('[OX2H,OX2R]')  # Include both free and bonded OH
    if not mol.HasSubstructMatch(hydroxyl):
        return False, "Missing required hydroxyl group"

    # Check for long carbon chains
    long_chain = Chem.MolFromSmarts('CCCCCCCC')  # At least 8 carbons in a chain
    if not mol.HasSubstructMatch(long_chain):
        return False, "Missing required long carbon chain"

    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 600:
        return False, "Molecular weight too low for sphingomyelin"

    # Check for sphingoid base connectivity
    # Pattern looks for -CH2-CH(NH)-CH(O)- backbone connected to phosphate
    sphingoid_base = Chem.MolFromSmarts('[CH2X4][CHX4]([NX3])[CHX4]([OX2,OX2H])COP(=O)')
    if not mol.HasSubstructMatch(sphingoid_base):
        # Try alternative pattern for different tautomers/representations
        alt_sphingoid = Chem.MolFromSmarts('[CH2X4][CX4]([NX3])[CX4]([OX2,OX2H])COP(=O)')
        if not mol.HasSubstructMatch(alt_sphingoid):
            return False, "Missing sphingoid base structure"

    return True, "Contains sphingoid base with amide-linked fatty acid and phosphocholine group"