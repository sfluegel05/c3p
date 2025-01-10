"""
Classifies: CHEBI:27300 vitamin D
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_vitamin_D(smiles: str):
    """
    Determines if a molecule is a vitamin D compound based on its SMILES string.
    Vitamin D compounds are defined as fat-soluble hydroxy seco-steroids,
    typically with a broken B-ring (seco-steroid) and specific sterochemistry.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a vitamin D compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more specific seco-steroid pattern with broken B-ring.
    # This pattern is adapted to capture the key structural features of vitamin D compounds.
    seco_steroid_pattern = Chem.MolFromSmarts("C1=CC[C@H]2CC[C@H](O)[C@H](C=3C[C@@H](O)CCC3=C)C2=C1") # Representative pattern focusing on seco-steroid

    if seco_steroid_pattern is None:
        return (None, None)  # Avoid failure if the pattern itself wasn't set up correctly

    if not mol.HasSubstructMatch(seco_steroid_pattern):
        return False, "Seco-steroid backbone typical of vitamin D not identified"

    # Improve the check for hydroxyl groups: vitamin D often has multiple hydroxylations
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]") # Detect OH groups
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 2:  # Expect at least two hydroxyl groups
        return False, f"Insufficient hydroxyl groups found, at least 2 expected"

    # Adjusted check to ensure the molecule is sufficiently hydrophobic/lipophilic
    mol_logP = rdMolDescriptors.CalcCrippenDescriptors(mol)[0]
    if mol_logP < 5:  # More typical threshold for lipophilicity in vitamin D
        return False, "Molecule not sufficiently hydrophobic for vitamin D classification"

    # Check molecular weight for vitamin D compounds
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_weight < 380:
        return False, "Molecular weight lower than typical vitamin D compounds"

    return True, "Matches expected vitamin D seco-steroid structure with sufficient hydroxyl groups and lipophilicity"