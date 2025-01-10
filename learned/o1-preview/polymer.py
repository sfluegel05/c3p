"""
Classifies: CHEBI:60027 polymer
"""
"""
Classifies: polymer
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdFMCS

def is_polymer(smiles: str):
    """
    Determines if a molecule is a polymer based on its SMILES string.
    A polymer is a macromolecule composed of repeating units connected in a chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polymer, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Remove counter ions and small fragments
    mol = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(mol) == 0:
        return False, "No valid molecular fragments found"
    elif len(mol) > 1:
        # Choose the largest fragment
        mol = max(mol, key=lambda m: m.GetNumAtoms())
    else:
        mol = mol[0]

    # Analyze the molecule for repeating units
    # Break single bonds to try to find repeating units
    bonds = mol.GetBonds()
    repeating_unit = None
    max_repeat_count = 0

    for bond in bonds:
        # Clone the molecule
        mol_clone = Chem.RWMol(mol)
        idx = bond.GetIdx()
        # Remove the bond
        mol_clone.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        # Get the fragments
        frags = Chem.GetMolFrags(mol_clone, asMols=True)
        if len(frags) < 2:
            continue
        # Try to find identical fragments
        frag_smiles = [Chem.MolToSmiles(frag) for frag in frags]
        frag_counts = {}
        for fs in frag_smiles:
            frag_counts[fs] = frag_counts.get(fs, 0) + 1
        # Identify the most common fragment
        common_frag = max(frag_counts, key=frag_counts.get)
        count = frag_counts[common_frag]
        if count > max_repeat_count:
            max_repeat_count = count
            repeating_unit = common_frag

    if repeating_unit and max_repeat_count >= 5:
        return True, f"Molecule has a repeating unit ({repeating_unit}) occurring {max_repeat_count} times, indicative of a polymer"

    return False, "No sufficient repeating units detected in the molecule"