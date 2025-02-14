"""
Classifies: CHEBI:16389 ubiquinones
"""
"""
Classifies: CHEBI:27022 ubiquinone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_ubiquinone(smiles: str):
    """
    Determines if a molecule is an ubiquinone based on its SMILES string.
    An ubiquinone is a benzoquinone derived from 2,3-dimethoxy-5-methylbenzoquinone,
    with a polyprenoid side chain typically attached at position 6.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ubiquinone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for 2,3-dimethoxy-5-methylbenzoquinone core
    core_pattern = Chem.MolFromSmarts("C1=C(C(=O)C(=C(C1=O)OC)OC)C")
    core_matches = mol.GetSubstructMatches(core_pattern)
    if not core_matches:
        return False, "Missing 2,3-dimethoxy-5-methylbenzoquinone core"

    # Look for polyprenoid side chain attached at position 6
    side_chain_pattern = Chem.MolFromSmarts("CC=C(C)CCC")
    side_chain_matches = mol.GetSubstructMatches(side_chain_pattern)
    if not side_chain_matches:
        return False, "No polyprenoid side chain found"

    # Check if side chain is attached to core at position 6
    for core_match in core_matches:
        for side_chain_match in side_chain_matches:
            if mol.GetBondBetweenAtoms(core_match[5], side_chain_match[0]).GetBondType() == Chem.BondType.SINGLE:
                return True, "Contains 2,3-dimethoxy-5-methylbenzoquinone core with polyprenoid side chain attached at position 6"

    return False, "Polyprenoid side chain not attached to the core at position 6"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:27022',
        'name': 'ubiquinone',
        'definition': 'Any benzoquinone derived from 2,3-dimethoxy-5-methylbenzoquinone; one of a group of naturally occurring homologues. The redox-active quinoid moiety usually carries a polyprenoid side chain at position 6, the number of isoprenoid units in which is species-specific. Ubiquinones are involved in the control of mitochondrial electron transport, and are also potent anti-oxidants.',
        'parents': ['CHEBI:33567', 'CHEBI:23888', 'CHEBI:36854']
    },
    'config': {...}
}