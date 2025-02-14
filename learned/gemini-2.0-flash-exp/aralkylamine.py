"""
Classifies: CHEBI:18000 aralkylamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aralkylamine(smiles: str):
    """
    Determines if a molecule is an aralkylamine based on its SMILES string.
    An aralkylamine is an alkylamine where at least one of the alkyl substituents is an aromatic group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aralkylamine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for an amine (primary, secondary, or tertiary, but not quaternary)
    amine_pattern = Chem.MolFromSmarts("[N;H0,H1,H2]") # Nitrogen with 0, 1 or 2 H
    if not mol.HasSubstructMatch(amine_pattern):
        return False, "No amine group found"

    # 2. Check for at least one aromatic ring
    aromatic_ring_pattern = Chem.MolFromSmarts("[c]")
    if not mol.HasSubstructMatch(aromatic_ring_pattern):
          return False, "No aromatic ring found"

    # 3. Check for an alkyl chain connecting the amine and an aromatic ring.
    # [N]-[C]-[#6]~[c] means: an amine connected to at least one carbon, followed by at least another atom connected to an aromatic carbon.
    alkyl_bridge_pattern = Chem.MolFromSmarts("[N]-[C]~[#6]~[c]")


    if not mol.HasSubstructMatch(alkyl_bridge_pattern):
        return False, "No alkyl chain connecting amine and aromatic ring"

    return True, "Contains an amine group connected to an aromatic ring via an alkyl chain"