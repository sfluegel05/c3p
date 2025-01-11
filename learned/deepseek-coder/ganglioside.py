"""
Classifies: CHEBI:28892 ganglioside
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ganglioside(smiles: str):
    """
    Determines if a molecule is a ganglioside based on its SMILES string.
    A ganglioside is a glycosphingolipid (ceramide and oligosaccharide) with one or more sialic acids linked to the sugar chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ganglioside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for ceramide structure (sphingosine + fatty acid)
    # Sphingosine pattern: C=C-C-C(OH)-C(NH)-C(OH)-C
    sphingosine_pattern = Chem.MolFromSmarts("[CH2X4][CHX4]=[CHX4][CHX4][CHX4]([OHX2])[CHX4]([NHX3])[CHX4]([OHX2])[CH2X4]")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No sphingosine backbone found (ceramide missing)"

    # Check for fatty acid (long carbon chain attached to NH2 of sphingosine)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 3:
        return False, "Missing fatty acid chain in ceramide"

    # Check for oligosaccharide (multiple sugar units)
    # Look for at least 3 sugar units (e.g., glucose, galactose, etc.)
    sugar_pattern = Chem.MolFromSmarts("[CX4][OX2][CX4]")
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if len(sugar_matches) < 3:
        return False, "Insufficient sugar units in oligosaccharide"

    # Check for sialic acid (N-acetylneuraminic acid, Neu5Ac)
    # Sialic acid pattern: C=O, NH, and multiple OH groups
    sialic_acid_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3][CX4][CX4]([OHX2])[CX4]([OHX2])[CX4]([OHX2])[CX4]([OHX2])")
    sialic_acid_matches = mol.GetSubstructMatches(sialic_acid_pattern)
    if len(sialic_acid_matches) < 1:
        return False, "No sialic acid (Neu5Ac) found"

    # Check molecular weight (gangliosides are typically large molecules)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for ganglioside"

    # Count carbons, oxygens, and nitrogens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    if c_count < 30:
        return False, "Too few carbons for ganglioside"
    if o_count < 10:
        return False, "Too few oxygens for ganglioside"
    if n_count < 2:
        return False, "Too few nitrogens for ganglioside"

    return True, "Contains ceramide, oligosaccharide, and sialic acid (ganglioside structure)"