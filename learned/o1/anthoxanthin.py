"""
Classifies: CHEBI:192499 anthoxanthin
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_anthoxanthin(smiles: str):
    """
    Determines if a molecule is an anthoxanthin based on its SMILES string.
    Anthoxanthins are flavonoid pigments, including flavones and flavonols,
    often with hydroxyl, methoxy, and glycosyl substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anthoxanthin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define flavone core (2-phenylchromen-4-one)
    flavone_core_smarts = 'O=C1C=CC2=CC=CC=C2O1-c3ccccc3'
    
    # Define flavonol core (3-hydroxyflavone)
    flavonol_core_smarts = 'O=C1C=CC2=CC=CC=C2O1-c3ccccc3O'
    
    # Create mol objects from SMARTS
    flavone_core_mol = Chem.MolFromSmarts(flavone_core_smarts)
    flavonol_core_mol = Chem.MolFromSmarts(flavonol_core_smarts)

    # Check for substructure match with flavone core or flavonol core
    is_flavone = mol.HasSubstructMatch(flavone_core_mol)
    is_flavonol = mol.HasSubstructMatch(flavonol_core_mol)

    # If match found, check for substitutions
    if is_flavone or is_flavonol:
        # Allow for hydroxyl and methoxy substitutions on the aromatic rings
        # Define patterns for hydroxyl and methoxy groups attached to aromatic rings
        substituent_smarts = '[cH]O | [cH]OC'
        substituent_mol = Chem.MolFromSmarts(substituent_smarts)

        # Check for substituents
        if mol.HasSubstructMatch(substituent_mol):
            return True, "Molecule is a substituted flavone or flavonol (anthoxanthin)"

        # Check for glycosides (sugar moieties)
        glycoside_smarts = '[O;!R][C@H]1O[C@H]([C@H]([C@H](O)[C@@H]1O)O)CO'
        glycoside_mol = Chem.MolFromSmarts(glycoside_smarts)
        if mol.HasSubstructMatch(glycoside_mol):
            return True, "Molecule is a glycosylated flavone or flavonol (anthoxanthin)"

        # If no substituents found, still consider as anthoxanthin
        return True, "Molecule contains flavone or flavonol core structure"

    # If no match is found
    return False, "Molecule does not match anthoxanthin criteria"