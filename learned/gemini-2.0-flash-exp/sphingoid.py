"""
Classifies: CHEBI:35785 sphingoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sphingoid(smiles: str):
    """
    Determines if a molecule is a sphingoid based on its SMILES string.
    A sphingoid has a long hydrocarbon chain, an amino group (primary, secondary, tertiary, or protonated), 
    and a vicinal hydroxyl or carbonyl group. Additionally, it may have a
    phosphocholine group or be glycosylated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sphingoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbons
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 10: # Sphingoid chain should be at least 10 carbons long
        return False, f"Carbon chain length too short: {carbon_count} carbons"
    
    # Check for various amine/amide groups
    amino_pattern1 = Chem.MolFromSmarts("[NX3H2]")  # primary amine
    amino_pattern2 = Chem.MolFromSmarts("[NH3+]")  # protonated primary amine
    amino_pattern3 = Chem.MolFromSmarts("[NX3H][CX4]") # secondary amine
    amino_pattern4 = Chem.MolFromSmarts("[NX3]([CX4])([CX4])") #tertiary amine
    amino_pattern5 = Chem.MolFromSmarts("[NX3][CX3](=[OX1])") # amide
    if not (mol.HasSubstructMatch(amino_pattern1) or mol.HasSubstructMatch(amino_pattern2) or
            mol.HasSubstructMatch(amino_pattern3) or mol.HasSubstructMatch(amino_pattern4) or
            mol.HasSubstructMatch(amino_pattern5)):
        return False, "No amine or amide group found"
    
    #Core sphingoid pattern check (C-C(OH/carbonyl)-C(N)) or variations
    # Pattern 1: C-C(OH)-C(N)-C
    sphingoid_pattern1 = Chem.MolFromSmarts("[CX4][CX4]([OH1])[CX4]([NX3])[CX4]")
    # Pattern 2: C-C(=O)-C(N)-C
    sphingoid_pattern2 = Chem.MolFromSmarts("[CX4][CX3](=[OX1])[CX4]([NX3])[CX4]")
     # Pattern 3: C-C(OH)-C(N)
    sphingoid_pattern3 = Chem.MolFromSmarts("[CX4][CX4]([OH1])[CX4]([NX3])")
    # Pattern 4: C-C(=O)-C(N)
    sphingoid_pattern4 = Chem.MolFromSmarts("[CX4][CX3](=[OX1])[CX4]([NX3])")
    #Pattern 5: C=C-C(OH)-C(N)
    sphingoid_pattern5 = Chem.MolFromSmarts("[CX4]=[CX3][CX4]([OH1])[CX4]([NX3])")
    #Pattern 6: C=C-C(=O)-C(N)
    sphingoid_pattern6 = Chem.MolFromSmarts("[CX4]=[CX3][CX4](=[OX1])[CX4]([NX3])")
    # Pattern 7: C-C(OH)-C(N)-C-OH
    sphingoid_pattern7 = Chem.MolFromSmarts("[CX4][CX4]([OH1])[CX4]([NX3])[CX4][OH1]")
     # Pattern 8: C-C(OH)-C(N)-C=C
    sphingoid_pattern8 = Chem.MolFromSmarts("[CX4][CX4]([OH1])[CX4]([NX3])[CX3]=[CX4]")
    # Pattern 9: C-C(=O)-C(N)-C=C
    sphingoid_pattern9 = Chem.MolFromSmarts("[CX4][CX3](=[OX1])[CX4]([NX3])[CX3]=[CX4]")
    # Pattern 10: C-C(OH)-C(N)-C-C-OH
    sphingoid_pattern10 = Chem.MolFromSmarts("[CX4][CX4]([OH1])[CX4]([NX3])[CX4][CX4][OH1]")

    if not (mol.HasSubstructMatch(sphingoid_pattern1) or mol.HasSubstructMatch(sphingoid_pattern2) or
            mol.HasSubstructMatch(sphingoid_pattern3) or mol.HasSubstructMatch(sphingoid_pattern4) or
            mol.HasSubstructMatch(sphingoid_pattern5) or mol.HasSubstructMatch(sphingoid_pattern6) or
            mol.HasSubstructMatch(sphingoid_pattern7) or mol.HasSubstructMatch(sphingoid_pattern8) or
            mol.HasSubstructMatch(sphingoid_pattern9) or mol.HasSubstructMatch(sphingoid_pattern10) ):
        return False, "Did not find the necessary hydroxyl/carbonyl and amine/amide groups with correct pattern"


    # Optional checks for phosphocholine and glycosylation
    phosphocholine_pattern = Chem.MolFromSmarts("COP([O-])(=O)OCC[N+](C)(C)C")
    glycoside_pattern = Chem.MolFromSmarts("OC[C@H]1[C@H](O)[C@H](O)[C@@H]([C@@H](O)O1)O") #check for glucose specifically
    
    
    is_phosphocholine = mol.HasSubstructMatch(phosphocholine_pattern)
    is_glycosylated = mol.HasSubstructMatch(glycoside_pattern)

    reason = "Contains a long hydrocarbon chain, an amino/amide group, and a vicinal hydroxyl/carbonyl group"
    if is_phosphocholine:
        reason += ", and phosphocholine group"
    if is_glycosylated:
        reason += ", and is glycosylated"
    
    return True, reason