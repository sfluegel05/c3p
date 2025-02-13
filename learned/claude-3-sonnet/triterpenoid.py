"""
Classifies: CHEBI:36615 triterpenoid
"""
"""
Classifies: CHEBI:36901 Triterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS
from rdkit.Chem.Descriptors import MolWt

def is_triterpenoid(smiles: str):
    """
    Determines if a molecule is a triterpenoid based on its SMILES string.
    Triterpenoids are terpenoids derived from a triterpene, where the C30 skeleton
    of the parent triterpene may have been rearranged or modified by the removal of
    skeletal atoms (generally methyl groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles, sanitize=True, removeHs=False)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for molecular weight in the range typical for triterpenoids
    mol_wt = MolWt(mol)
    if mol_wt < 400 or mol_wt > 600:
        return False, f"Molecular weight {mol_wt:.2f} outside typical triterpenoid range"

    # Check for the presence of certain substructures/functional groups
    substructures = [
        Chem.MolFromSmarts("[C&ring]1[C&ring][C&ring][C&ring][C&ring][C&ring][C&ring]1"),  # Cyclic structures
        Chem.MolFromSmarts("[OH]"),  # Hydroxyl groups
        Chem.MolFromSmarts("[C&ring](=O)"),  # Cyclic carbonyl groups
        Chem.MolFromSmarts("[C&ring]1[C&ring][C&ring][C&ring][C&ring][C&ring][C&ring]1"),  # Additional cyclic structures
    ]
    
    for substructure in substructures:
        if not mol.HasSubstructMatch(substructure):
            return False, f"Missing substructure: {Chem.MolToSmarts(substructure)}"

    # Check for specific triterpene skeletons
    skeletons = [
        Chem.MolFromSmiles("C[C@H]1[C@H]2[C@@H](C[C@@H]3[C@@H]4CC[C@H]5C[C@@H](O)CC[C@]5(C)[C@H]4CC[C@]23C)O[C@@]11"),  # Oleanane
        Chem.MolFromSmiles("C[C@@H]1[C@@H]2[C@@H](C[C@@H]3[C@@H]4CC[C@H]5C[C@@H](O)CC[C@]5(C)[C@H]4CC[C@]23C)O[C@@]11"),  # Ursane
        Chem.MolFromSmiles("C[C@@H]1[C@H]2[C@@H](C[C@@H]3[C@@H]4CC[C@H]5C[C@@H](O)CC[C@]5(C)[C@H]4CC[C@]23C)O[C@@]11"),  # Lupane
        Chem.MolFromSmiles("C[C@@H]1[C@H]2[C@@H](C[C@@H]3[C@@H]4CC[C@H]5C[C@@H](O)CC[C@]5(C)[C@H]4CC[C@]23C)O[C@@]11"),  # Dammarane
    ]

    skeleton_match = False
    for skeleton in skeletons:
        mcs = rdFMCS.FindMCS([mol, skeleton], matchValences=True, completeRingsOnly=True)
        if mcs.numAtoms > 20:  # Adjust threshold as needed
            skeleton_match = True
            break

    if not skeleton_match:
        return False, "No known triterpene skeleton found"

    return True, "Molecule exhibits characteristics of a triterpenoid"