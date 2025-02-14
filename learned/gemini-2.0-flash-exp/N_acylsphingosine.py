"""
Classifies: CHEBI:52639 N-acylsphingosine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_N_acylsphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylsphingosine based on its SMILES string.
    An N-acylsphingosine has a sphingosine backbone with a fatty acyl group attached to the nitrogen via an amide bond.

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

    # Define SMARTS for sphingosine backbone (more generic, ignoring stereochemistry)
    # [CX4,CX3]-[CX4,CX3](O)-[CX4,CX3](N)-[CX4,CX3](O)-[CX4,CX3]=[CX4,CX3]
    sphingosine_pattern = Chem.MolFromSmarts("[CX4,CX3][CX4,CX3](O)[CX4,CX3](N)[CX4,CX3](O)[CX4,CX3]=[CX4,CX3]")
    if not mol.HasSubstructMatch(sphingosine_pattern):
       return False, "No sphingosine backbone found"

    # Define SMARTS for the amide bond
    amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)

    if len(amide_matches) !=1:
        return False, f"Found {len(amide_matches)} amide bonds, need exactly 1"

    # Check for long carbon chain by checking the number of carbons
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 18:
       return False, "Too few carbons for N-acylsphingosine"

    # Check for at least one hydroxyl group
    oh_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and len(atom.GetNeighbors()) == 1)
    if oh_count < 1:
        return False, "Must have at least one hydroxyl group"
    
    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
      return False, "Molecular weight too low for N-acylsphingosine"
    
    return True, "Contains a sphingosine backbone with a fatty acyl group attached via an amide bond"