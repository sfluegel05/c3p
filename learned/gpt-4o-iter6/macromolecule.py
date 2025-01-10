"""
Classifies: CHEBI:33839 macromolecule
"""
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_macromolecule(smiles: str):
    """
    Determines if a molecule is a macromolecule based on its SMILES string.
    A macromolecule is defined as a molecule with high molecular mass,
    composed of repetitive units derived from low molecular mass molecules.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a macromolecule, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    # Check for high molecular weight (arbitrary cutoff for demonstration, e.g., 2000 Da)
    if mol_wt >= 2000:
        return True, f"Molecule has high molecular weight: {mol_wt:.2f} Da"
    
    # Approximate check for repeating units by counting large ring systems or long chains
    # Check number of heavy atoms
    num_atoms = mol.GetNumAtoms()
    
    # Some polymers might be linear with a high number of heavy atoms
    if num_atoms >= 100:  # Arbitrary number for demonstration
        return True, f"Molecule has a large number of atoms: {num_atoms}"

    return False, "Molecule does not meet criteria for a macromolecule"

# Example testing
example_smiles = "C1CC2=C(NC1=O)C(=O)N(C(=O)N2)C(C(=O)N3CCC(=O)N(C3=O)C4=CC=CC=C4)C5CCCN5"
print(is_macromolecule(example_smiles))