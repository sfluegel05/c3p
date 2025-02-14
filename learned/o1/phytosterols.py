"""
Classifies: CHEBI:26125 phytosterols
"""
"""
Classifies: CHEBI:18374 phytosterol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phytosterols(smiles: str):
    """
    Determines if a molecule is a phytosterol based on its SMILES string.
    A phytosterol is a plant sterol similar to cholesterol, varying only in carbon side chains
    and/or presence or absence of double bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phytosterol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Steroid nucleus SMILES (cyclopentanoperhydrophenanthrene skeleton)
    sterol_core_smiles = "C1CCC2C(C1)CCC3C2CCC4C3CCCC4"
    sterol_core = Chem.MolFromSmiles(sterol_core_smiles)

    # Check if molecule contains sterol core
    if not mol.HasSubstructMatch(sterol_core):
        return False, "Steroid nucleus not found"

    # Optional: Exclude cholesterol by checking the side chain at position 17
    # Phytosterols have variations in the side chain compared to cholesterol
    # We can check the substituent at C17 (the carbon attached to ring D)

    # Get the atom index of C17 in the sterol core
    match = mol.GetSubstructMatch(sterol_core)
    if not match:
        return False, "Steroid nucleus not found"
    c17_index = match[16]  # Indexing starts from 0

    # Get the substituents attached to C17
    c17_atom = mol.GetAtomWithIdx(c17_index)
    neighbors = c17_atom.GetNeighbors()
    side_chain_atoms = [atom for atom in neighbors if atom.GetIdx() not in match]

    # Check if there is a side chain attached to C17
    if not side_chain_atoms:
        return False, "No side chain attached at position C17"

    # Optionally, check the size of the side chain (phytosterols typically have longer side chains)
    side_chain_atom = side_chain_atoms[0]  # Assuming one side chain
    side_chain_frag = Chem.PathToSubmol(mol, Chem.rdmolops.GetShortestPath(mol, c17_index, side_chain_atom.GetIdx()))
    side_chain_mw = Chem.Descriptors.ExactMolWt(side_chain_frag)

    # If needed, set a threshold for side chain molecular weight
    if side_chain_mw < 40:
        return False, "Side chain at C17 is too small for phytosterol"

    return True, "Molecule contains steroid nucleus characteristic of phytosterols"