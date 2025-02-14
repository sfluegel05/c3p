"""
Classifies: CHEBI:26199 polyprenol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_polyprenol(smiles: str):
    """
    Determines if a molecule is a polyprenol based on its SMILES string.
    Polyprenols are characterized by a repeating isoprene units, arranged in a
    linear chain, terminated by an alcohol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyprenol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for a single isoprene unit within a polyprenol chain, including the methyl group and single bonds between units.
    # The pattern looks for a head-to-tail connection of isoprene units -C=C(C)-C-C-. It assumes the double bond is in the E configuration.
    # The asterisk (*) matches any atom, ensuring it can be part of a chain. The double bond has to be connected to a carbon, and that carbon must be connected to a CH3
    isoprene_unit_pattern = Chem.MolFromSmarts("[*]C(=C([CH3])[CH2])[CH2][*]")


    matches = mol.GetSubstructMatches(isoprene_unit_pattern)
    num_isoprene_units = len(matches)

    if num_isoprene_units < 2:
      return False, f"Too few isoprene units, found {num_isoprene_units}"

    # Check for terminal alcohol at the end of chain (using SMARTS, which also checks connection to C=C)
    terminal_alcohol_pattern = Chem.MolFromSmarts("[CH2][CH2][C](=[CH])[CH2][OX2H]") #terminal alcohol has to be connected to an isoprene unit
    terminal_alcohol_matches = mol.GetSubstructMatches(terminal_alcohol_pattern)

    if len(terminal_alcohol_matches) != 1:
        return False, "Incorrect number of terminal alcohol groups. Should be exactly 1"

    # Additional check: Ensure there are not any other alcohols within the chain
    all_alcohol_pattern = Chem.MolFromSmarts("[CX4][OX2H]") # all alcohols
    all_alcohol_matches = mol.GetSubstructMatches(all_alcohol_pattern)

    if len(all_alcohol_matches) > 1:
        return False, "More than one alcohol group detected, which is not compatible with polyprenols"


    # Check for chain length (number of carbon atoms) - must be at least 4*num_isoprene_units, + one for the terminal alcohol
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < (4 * num_isoprene_units + 1):
      return False, "Too few carbons for the number of isoprene units."

    # Check rotatable bonds (expect at least 2 * num_isoprene_units, - 1 cause not all single bonds are rotatable)
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < (2 * num_isoprene_units-1) :
        return False, "Too few rotatable bonds for polyprenol chain."


    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 150: # minimum mass of 2 isoprene and one alcohol
      return False, "Molecular weight is too low"

    return True, "Molecule has repeating isoprene units, a terminal alcohol, sufficient chain length, and is compatible with polyprenol definition."