"""
Classifies: CHEBI:33839 macromolecule
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_macromolecule(smiles: str):
    """
    Determines if a molecule is a macromolecule based on its SMILES string.
    Due to the diversity and complexity of macromolecules, the classification
    relies on general features such as large molecular size and complex structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule could be considered a macromolecule, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Calculate the number of atoms and molecular weight
    num_atoms = mol.GetNumAtoms()
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    
    # Using a general molecular weight threshold for macromolecules
    if num_atoms > 50 and mol_weight > 1000:
        return True, "Large number of atoms and high molecular weight indicate a macromolecule"
    
    return False, "Does not meet general macromolecule criteria"