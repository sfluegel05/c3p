"""
Classifies: CHEBI:33839 macromolecule
"""
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_macromolecule(smiles: str):
    """
    Determines if a molecule is a macromolecule based on its SMILES string.
    A macromolecule is a molecule of high relative molecular mass, with a structure
    that essentially comprises the multiple repetition of units derived from smaller molecules.

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
    mol_wt = Descriptors.MolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight ({mol_wt:.2f} Da) is too low for a macromolecule"

    # Check for repeating units (e.g., polymers)
    # This is a simplistic check; more sophisticated methods could be used
    repeating_pattern = Chem.MolFromSmarts("[*]~[*]~[*]~[*]")
    matches = mol.GetSubstructMatches(repeating_pattern)
    if len(matches) < 3:
        return False, "No significant repeating units found"

    # Check for complexity (number of atoms and bonds)
    num_atoms = mol.GetNumAtoms()
    num_bonds = mol.GetNumBonds()
    if num_atoms < 50 or num_bonds < 50:
        return False, f"Not complex enough (atoms: {num_atoms}, bonds: {num_bonds})"

    # Check for high number of rotatable bonds (indicative of flexibility in macromolecules)
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, f"Too few rotatable bonds ({n_rotatable}) for a macromolecule"

    return True, "High molecular weight, repeating units, and complexity suggest a macromolecule"