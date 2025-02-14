"""
Classifies: CHEBI:31998 N-acylphytosphingosine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_N_acylphytosphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylphytosphingosine based on its SMILES string.
    An N-acylphytosphingosine is a phytosphingosine with a fatty acyl group attached to the nitrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylphytosphingosine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for phytosphingosine backbone using SMARTS
    #   - The pattern checks for C-C(O)-C(O)-C(N)-C(O) where the N atom is connected to three other atoms (to identify the N in the chain)
    #   - The 3 chiral centers are not enforced. This avoids errors in molecules that have stereoisomers or unspecified stereocenters.
    phytosphingosine_pattern = Chem.MolFromSmarts("[CH2X4][CHX4]([OX2])[CHX4]([OX2])[CHX4]([NX3])[CHX4]([OX2])") 
    if not mol.HasSubstructMatch(phytosphingosine_pattern):
       return False, "Phytosphingosine backbone not found"

    # 2. Check for amide linkage (-C(=O)N-)
    amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3]")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) < 1:
        return False, "No amide linkage found"
    
    #3. Check for fatty acyl chain
    #   -  Looking for a long carbon chain connected to the carbonyl of the amide linkage
    fatty_acyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") # minimum chain of 4 C atoms
    fatty_acyl_matches = mol.GetSubstructMatches(fatty_acyl_pattern)
    
    if len(fatty_acyl_matches) < 1:
      return False, "Fatty acyl chain not found"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Acyl chain too short to be a fatty acid"

    # Check molecular weight - N-acylphytosphingosines typically > 400 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for N-acylphytosphingosine"
    
    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)

    if c_count < 18:
        return False, "Too few carbons for N-acylphytosphingosine"
    if o_count < 4:
        return False, "Too few oxygens for N-acylphytosphingosine"
    if n_count != 1:
      return False, "Must have one nitrogen in the molecule"

    return True, "Contains a phytosphingosine backbone with a fatty acyl group attached to the nitrogen"