"""
Classifies: CHEBI:72588 semisynthetic derivative
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Descriptors

def is_semisynthetic_derivative(smiles: str):
    """
    Determines if a molecule could be classified as a semisynthetic derivative based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely a semisynthetic derivative, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Look for typical semisynthetic modification groups
    # Ester groups in general
    ester_group = Chem.MolFromSmarts("C(=O)O")
    if mol.HasSubstructMatch(ester_group):
        return True, "Contains ester group, indicates possible chemical modification"

    # Ether bonds
    ether_group = Chem.MolFromSmarts("COC")
    if mol.HasSubstructMatch(ether_group):
        return True, "Contains ether bond, indicates possible chemical modification"

    # Tertiary amines (methylated amines)
    tertiary_amine = Chem.MolFromSmarts("N(C)(C)C")
    if mol.HasSubstructMatch(tertiary_amine):
        return True, "Contains tertiary amine, common in semisynthetic derivatives"
    
    # Chlorination/Bromination (halogenated aromatics)
    halogenated_aromatic = Chem.MolFromSmarts("c[F,Cl,Br,I]")
    if mol.HasSubstructMatch(halogenated_aromatic):
        return True, "Contains halogenated aromatic, indicates possible chemical modification"

    # Amide groups extensively used in derivatization
    amide_group = Chem.MolFromSmarts("C(=O)N")
    if mol.HasSubstructMatch(amide_group):
        return True, "Contains amide group, common in semisynthetic derivatives"

    # Look for branching (indicating complex modifications)
    branching_pattern = Chem.MolFromSmarts("C(C)(C)C")
    if mol.HasSubstructMatch(branching_pattern):
        return True, "Contains branching pattern, indicative of derivatization"

    return False, "Molecule does not have typical indicators of semisynthetic derivation"