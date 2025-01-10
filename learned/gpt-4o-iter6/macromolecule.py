"""
Classifies: CHEBI:33839 macromolecule
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_macromolecule(smiles: str):
    """
    Determines if a molecule is a macromolecule based on its SMILES string.
    A macromolecule is defined as a molecule with high molecular mass,
    composed of repetitive units derived from low molecular mass molecules.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a macromolecule, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    # Check for high molecular weight - typically in 5000 Da or greater for common macromolecules
    if mol_wt >= 5000:
        return True, f"Molecule has high molecular weight: {mol_wt:.2f} Da"
    
    # Count total number of heavy atoms (typically exceeding 100 indicates macromolecule)
    atom_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)
    
    # For checking repeating units, we look at ring or similar structures — indicative for synthetic polymers like nylons, polyesters
    unique_subunits = Chem.GetSymmSSSR(mol)
    
    # Approximate check for large number of atoms and repeating subunits might denote a polymeric nature
    if atom_count >= 100 or len(unique_subunits) > 10:
        return True, f"Molecule has a large number of atoms ({atom_count}) or repeating substructures ({len(unique_subunits)})"
    
    return False, "Molecule does not meet criteria for a macromolecule"

# Example testing
example_smiles = "O1[C@H](O[C@H]2OC(CO)C(O)C2O[C@@H](C(O)CO)O[C@@H]3O[C@@H]([C@@H](O)C(O)[C@H]3O)CO)C(O)C(O)[C@H](CO)O[C@H]4OC(CO)C(O)C4O"
print(is_macromolecule(example_smiles))