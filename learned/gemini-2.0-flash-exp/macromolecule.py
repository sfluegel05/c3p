"""
Classifies: CHEBI:33839 macromolecule
"""
"""
Classifies: CHEBI:25367 macromolecule
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_macromolecule(smiles: str):
    """
    Determines if a molecule is a macromolecule based on its SMILES string.
    A macromolecule is a large molecule built from repeating smaller units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a macromolecule, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)

    # Check number of heavy (non-H) atoms
    num_heavy_atoms = mol.GetNumHeavyAtoms()

    # Calculate number of rotatable bonds - indicates chain length
    num_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    
    # Apply heuristics for macromolecule classification
    if mol_wt > 800 and num_heavy_atoms > 40 and num_rotatable_bonds > 10:
          return True, "High molecular weight, large number of heavy atoms and rotatable bonds indicating a macromolecule."
    elif mol_wt > 1200:
          return True, "Very high molecular weight, likely a macromolecule."
    elif num_heavy_atoms > 60 and num_rotatable_bonds > 20:
          return True, "Large number of heavy atoms and rotatable bonds, indicating a macromolecule"
    else:
        return False, "Does not meet macromolecule criteria based on molecular weight, heavy atoms and number of rotatable bonds."