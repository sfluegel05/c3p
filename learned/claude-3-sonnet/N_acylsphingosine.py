"""
Classifies: CHEBI:52639 N-acylsphingosine
"""
"""
Classifies: CHEBI:37713 N-acylsphingosine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_N_acylsphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylsphingosine based on its SMILES string.
    An N-acylsphingosine contains a sphingosine backbone with an unspecified fatty acyl group attached to the nitrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylsphingosine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for sphingosine backbone pattern
    sphingosine_pattern = Chem.MolFromSmarts("[C@@H]([C@@H](/C=C/CCCCCCC)O)CO")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No sphingosine backbone found"

    # Look for amide group (-N-C(=O)-) attached to sphingosine backbone
    amide_pattern = Chem.MolFromSmarts("[N]([CH2][CH]([OH])[CH](/C=C/CCCCCCC))[C](=O)")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) == 0:
        return False, "No amide group attached to sphingosine backbone"

    # Look for fatty acid chain attached to amide group
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern, maxMatches=1, useChirality=True)
    if len(fatty_acid_matches) == 0:
        return False, "No fatty acid chain attached to amide group"

    # Check for additional functional groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    
    linoleoyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2]CCCCCCC/C=C\C/C=C\CCCCC")
    linoleoyl_matches = mol.GetSubstructMatches(linoleoyl_pattern)

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Chains too short to be fatty acids"

    # Adjust molecular weight and atom count filters
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400 or mol_wt > 1000:
        return False, "Molecular weight outside expected range for N-acylsphingosine"

    # Check for reasonable atom counts
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    if c_count < 20 or o_count < 3 or n_count != 1:
        return False, "Atom counts outside expected range for N-acylsphingosine"

    # Construct reason string based on additional functional groups
    reason = "Contains sphingosine backbone with fatty acyl chain attached to nitrogen"
    if len(hydroxyl_matches) > 0:
        reason += ", with additional hydroxyl group(s)"
    if len(linoleoyl_matches) > 0:
        reason += ", with linoleoyl group"

    return True, reason