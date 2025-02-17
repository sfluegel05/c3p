"""
Classifies: CHEBI:33916 aldopentose
"""
"""
Classifies: aldopentose, defined as 'A pentose with a (potential) aldehyde group at one end.'
Examples include various forms of xylose, ribose, arabinose, lyxose, etc.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_aldopentose(smiles: str):
    """
    Determines if a molecule is an aldopentose based on its SMILES string.
    An aldopentose is a five‐carbon sugar that either shows an aldehyde group in the open-chain form 
    or exists in a cyclic form (hemiacetal) capable of generating an aldehyde upon ring opening.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is an aldopentose, False otherwise.
        str: The reason for the classification.
    """
    
    # Parse SMILES string to get a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count the number of carbon atoms in the molecule.
    # An aldopentose should have exactly 5 carbon atoms.
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(carbon_atoms) != 5:
        return False, f"Number of carbon atoms is {len(carbon_atoms)}; expected exactly 5 for a pentose"
    
    # Count the number of oxygen atoms.
    # For a typical aldopentose (C5H10O5), there should be exactly 5 oxygen atoms.
    oxygen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8]
    if len(oxygen_atoms) != 5:
        return False, f"Number of oxygen atoms is {len(oxygen_atoms)}; expected 5 for C5H10O5"
    
    # Count hydroxyl groups (-OH) using a SMARTS pattern.
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")  # oxygen with one hydrogen
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 3:
        return False, f"Found only {len(hydroxyl_matches)} -OH groups; expected at least 3 for an aldopentose"
    
    # Check for an aldehyde group in the open-chain form.
    # Aldehyde group SMARTS: a carbon with one hydrogen bound to =O.
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)")
    if mol.HasSubstructMatch(aldehyde_pattern):
        return True, "Open-chain aldopentose: aldehyde group detected."
    else:
        # If an aldehyde group is not detected, check if the molecule is cyclic.
        # In cyclic forms the aldehyde group is 'hidden' as a hemiacetal.
        if mol.GetRingInfo().NumRings() > 0:
            return True, "Cyclized aldopentose: potential open-chain aldehyde form upon ring opening."
        else:
            return False, "No aldehyde group detected and the molecule is not cyclic; does not match aldopentose criteria."
    
# (Optional: for testing, one might call the function with one of the sample SMILES)
# e.g., print(is_aldopentose("[H]C(=O)[C@@H](O)[C@H](O)[C@@H](O)CO"))