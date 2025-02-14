"""
Classifies: CHEBI:26267 proanthocyanidin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_proanthocyanidin(smiles: str):
    """
    Determines if a molecule is a proanthocyanidin based on its SMILES string.
    A proanthocyanidin is a flavonoid oligomer formed by condensation of two or more hydroxyflavan units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a proanthocyanidin, False otherwise.
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a general flavan unit pattern
    flavan_pattern = Chem.MolFromSmarts("[c]1[c]([OH])[c]([OH])[c](C[C]([OH])C2[O][c]3[c]([OH])[c]([OH])[c][c]2[c]3)[c][c]1")
    if flavan_pattern is None:
        return False, "Invalid flavan unit SMARTS pattern"

    flavan_matches = mol.GetSubstructMatches(flavan_pattern)
    num_flavan_units = len(flavan_matches)

    if num_flavan_units < 2:
        return False, f"Found {num_flavan_units} flavan units, need at least 2."

    # Additional check: Count the number of aromatic rings. Proanthocyanidins should have a minimum of 3
    aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if aromatic_rings < 3:
        return False, "Too few aromatic rings for a proanthocyanidin"

        # Check the molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:  # Lowered the threshold
      return False, "Too small molecular weight for a proanthocyanidin"

    # Count carbons and oxygens to verify
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 20:
        return False, "Too few carbons for proanthocyanidin"
    if o_count < 8:
        return False, "Too few oxygens for proanthocyanidin"

    return True, "Contains at least two linked flavan units."