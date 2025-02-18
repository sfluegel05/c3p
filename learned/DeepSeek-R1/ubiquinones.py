"""
Classifies: CHEBI:16389 ubiquinones
"""
"""
Classifies: CHEBI:46245 ubiquinones
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_ubiquinones(smiles: str):
    """
    Determines if a molecule is a ubiquinone based on its SMILES string.
    Ubiquinones are derived from 2,3-dimethoxy-5-methylbenzoquinone and typically have a polyprenoid side chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ubiquinone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define core benzoquinone pattern: 2,3-dimethoxy-5-substituted-1,4-benzoquinone
    # SMARTS pattern for the core structure (positions may vary based on SMILES numbering)
    core_pattern = Chem.MolFromSmarts("[O]=C1C(OC)=C(OC)C(=O)C(-[#6])C=C1")
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Core 2,3-dimethoxy-5-substituted benzoquinone structure not found"

    # Optional check for polyprenoid side chain (long isoprenoid chain)
    # Look for at least 3 consecutive isoprene units (C=C-C-C)
    isoprene_pattern = Chem.MolFromSmarts("C=CC(CC=CC(CC=CCCC)CC)")
    if mol.HasSubstructMatch(isoprene_pattern):
        return True, "Core benzoquinone with polyprenoid side chain"
    
    # Accept molecules without long chains (e.g. ubiquinone-0 has methyl group)
    return True, "Core benzoquinone structure present (may have short side chain)"