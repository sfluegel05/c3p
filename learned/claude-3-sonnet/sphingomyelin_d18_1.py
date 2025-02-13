"""
Classifies: CHEBI:17636 sphingomyelin d18:1
"""
"""
Classifies: sphingomyelin d18:1
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sphingomyelin_d18_1(smiles: str):
    """
    Determines if a molecule is a sphingomyelin d18:1 based on its SMILES string.
    Sphingomyelin d18:1 has a sphingosine backbone (18 carbons, 1 double bond),
    with specific stereochemistry and a phosphocholine head group.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a sphingomyelin d18:1, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phosphocholine head group
    phosphocholine = Chem.MolFromSmarts("[O-]P(=O)(OCC[N+](C)(C)C)O")
    if not mol.HasSubstructMatch(phosphocholine):
        return False, "Missing phosphocholine head group"

    # Check for sphingosine backbone with correct stereochemistry
    # [C@H] ensures R stereochemistry at the 2-position (amino)
    # [C@@H] ensures R stereochemistry at the 3-position (hydroxyl)
    # \C=C\ ensures trans/E geometry of the double bond
    sphingosine = Chem.MolFromSmarts("[CH2]OP([O-])(=O)OCC[N+](C)(C)C.[C@H](CN)([C@@H](O)/C=C/CCCCCCCCCCCCC)")
    if not mol.HasSubstructMatch(sphingosine):
        return False, "Missing or incorrect sphingosine d18:1 backbone"

    # Check for amide linkage
    amide = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[#6]")
    amide_matches = mol.GetSubstructMatches(amide)
    if len(amide_matches) != 1:
        return False, "Must have exactly one amide bond"

    # Verify presence of long chain fatty acid
    fatty_acid = Chem.MolFromSmarts("C(=O)CCCCCCCC")  # At least 8 carbons in acyl chain
    if not mol.HasSubstructMatch(fatty_acid):
        return False, "Missing long chain fatty acid"

    # Count key atoms
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if n_count != 2:  # One from amide, one from choline
        return False, "Must have exactly 2 nitrogen atoms"
    if p_count != 1:
        return False, "Must have exactly 1 phosphorus atom"
    if o_count < 6:  # Minimum 6 oxygens (phosphate, hydroxyl, amide, phosphocholine)
        return False, "Insufficient oxygen atoms"

    # Verify molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 600:
        return False, "Molecular weight too low for sphingomyelin"

    # Count carbons to verify total size
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 35:  # Minimum expected carbons (18 from sphingosine + ~14 from fatty acid + 5 from phosphocholine)
        return False, "Insufficient carbon atoms for sphingomyelin"

    # Additional check for trans double bond in sphingosine
    trans_alkene = Chem.MolFromSmarts("C/C=C/C")
    if not mol.HasSubstructMatch(trans_alkene):
        return False, "Missing required trans double bond"

    return True, "Contains sphingosine d18:1 backbone with correct stereochemistry, phosphocholine head group, and N-acyl chain"