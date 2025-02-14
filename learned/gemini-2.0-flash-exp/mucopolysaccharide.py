"""
Classifies: CHEBI:37395 mucopolysaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_mucopolysaccharide(smiles: str):
    """
    Determines if a molecule is likely a mucopolysaccharide based on its SMILES string using heuristics.
    
    Refined version with more specific substructure checks and molecular weight considerations.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple: (bool, str) - True if likely a mucopolysaccharide, False otherwise, along with a reason.
               Returns (None, None) if SMILES processing fails.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None # Indicate processing failure

    # Define SMARTS patterns
    # Uronic acid (a six-membered ring with a carboxyl group at C6)
    uronic_acid_pattern = Chem.MolFromSmarts("[CX1H2]([OX2H])[CX4H]([OX2H])[CX4H]([OX2H])[CX4H]([OX2H])[CX4H]([OX2H])[CX3](=O)[OX1H]")
    uronic_acid_matches = mol.GetSubstructMatches(uronic_acid_pattern)
    
    # Sulfate on a sugar ring
    sulfate_sugar_pattern1 = Chem.MolFromSmarts("[OX2H][CX4]([OX2H])([OX2H])[CX4]([OX2H])([OX2H])[CX4]([OX2H])([OX2H])[CX4]([OX2H])([OX2H])[CX4]([OX2H])([OX2H])[OX2]S(=O)(=O)O")
    sulfate_sugar_pattern2 = Chem.MolFromSmarts("[CX4][OX2]S(=O)(=O)O") #Sulfate on ring
    sulfate_matches1 = mol.GetSubstructMatches(sulfate_sugar_pattern1)
    sulfate_matches2 = mol.GetSubstructMatches(sulfate_sugar_pattern2)
    sulfate_matches = sulfate_matches1 + sulfate_matches2

    # Glycosamine (amino group on a sugar ring)
    glycosamine_pattern = Chem.MolFromSmarts("[CX1H2]([OX2H])[CX4H]([OX2H])[CX4H]([OX2H])[CX4H]([OX2H])[CX4H]([OX2H])[NX3]")
    glycosamine_matches = mol.GetSubstructMatches(glycosamine_pattern)
    
    # Basic amino group (for a general check as aminosugar may be a bad SMARTS)
    amino_pattern = Chem.MolFromSmarts("[NX3]")
    amino_matches = mol.GetSubstructMatches(amino_pattern)

    # Basic oxygen (only ring or non-carbonyl oxygens)
    oxygen_pattern = Chem.MolFromSmarts("[OX2]")
    oxygen_matches = mol.GetSubstructMatches(oxygen_pattern)

    carbon_pattern = Chem.MolFromSmarts("[CX4]")
    carbon_matches = mol.GetSubstructMatches(carbon_pattern)
    
    n_atoms = mol.GetNumAtoms()
    if n_atoms < 10:
        return False, "Too small for mucopolysaccharide"

    # Check for minimal number of uronic/sulfate and glycosamine units.
    if len(uronic_acid_matches) < 1 and len(sulfate_matches) < 1:
        return False, "No uronic acid or sulfate groups found."
    if len(glycosamine_matches) < 1 and len(amino_matches) <1:
        return False, "No glycosamine or amino groups found."

    #Check for a good ratio of oxygens to carbons, avoid small, highly functionalized molecules
    if len(oxygen_matches) < 5 :
      return False, "Too few oxygen atoms"
    if len(carbon_matches) < 10:
      return False, "Too few carbons"
    if len(oxygen_matches) > 2*len(carbon_matches):
        return False, "Too many oxygen atoms"

    # Check molecular weight and number of rotatable bonds.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for a mucopolysaccharide."

    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Too few rotatable bonds for a mucopolysaccharide."
    
    return True, "Likely mucopolysaccharide: contains uronic/sulfate groups, amino groups, and has sufficient mass/rotatable bonds."