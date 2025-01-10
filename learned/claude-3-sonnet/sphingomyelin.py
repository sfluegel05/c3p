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

    # Check for phosphocholine group with more flexible pattern
    # Allow for different representations of the phosphocholine group
    phosphocholine_patterns = [
        '[OX2]P(=O)([OX2-])OCC[NX4+](C)(C)C',
        '[OX2]P([OX2-])(=O)OCC[NX4+](C)(C)C'
    ]
    has_phosphocholine = False
    for pattern in phosphocholine_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_phosphocholine = True
            break
    if not has_phosphocholine:
        return False, "Missing phosphocholine group"

    # Check for amide linkage (-NC(=O)-)
    amide = Chem.MolFromSmarts('[NX3][CX3](=[OX1])[#6]')
    if not mol.HasSubstructMatch(amide):
        return False, "Missing amide linkage"

    # Check for sphingoid base characteristics
    # Look for either sphingosine (with double bond) or sphinganine (without double bond) backbone
    sphingoid_patterns = [
        # Pattern for sphingosine-type backbone (with double bond)
        '[OX2H][CX4][CX4]([NX3])[CX4]O[PX4]',
        # Pattern for sphinganine-type backbone (saturated)
        '[OX2H][CX4][CX4]([NX3])[CX4]O[PX4]',
    ]
    has_sphingoid = False
    for pattern in sphingoid_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_sphingoid = True
            break
    if not has_sphingoid:
        return False, "Missing sphingoid base structure"

    # Count key atoms to verify composition
    atom_counts = {}
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        atom_counts[symbol] = atom_counts.get(symbol, 0) + 1

    # Verify minimum requirements for sphingomyelin
    if atom_counts.get('C', 0) < 20:
        return False, "Insufficient carbon atoms for sphingomyelin"
    if atom_counts.get('N', 0) != 2:
        return False, "Must have exactly 2 nitrogen atoms"
    if atom_counts.get('P', 0) != 1:
        return False, "Must have exactly 1 phosphorus atom"
    if atom_counts.get('O', 0) < 4:
        return False, "Insufficient oxygen atoms"

    # Check for long carbon chains (characteristic of sphingomyelin)
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Insufficient chain length for sphingomyelin"

    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:  # Lowered from 600 to be more inclusive
        return False, "Molecular weight too low for sphingomyelin"

    # Additional check for hydroxyl group
    hydroxyl = Chem.MolFromSmarts('[OX2H]')
    if not mol.HasSubstructMatch(hydroxyl):
        return False, "Missing required hydroxyl group"

    return True, "Contains sphingoid base with amide-linked fatty acid and phosphocholine group"