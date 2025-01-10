"""
Classifies: CHEBI:18310 alkane
"""
"""
Classifies: alkane
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alkane(smiles: str):
    """
    Determines if a molecule is an alkane based on its SMILES string.
    An alkane is an acyclic branched or unbranched hydrocarbon having the general formula CnH2n+2,
    consisting entirely of hydrogen atoms and saturated carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkane, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for atoms other than carbon and hydrogen
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num != 6 and atomic_num != 1:
            return False, f"Contains atom other than C and H: {atom.GetSymbol()}"

    # Check for rings (molecule must be acyclic)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Contains rings, molecule is cyclic"

    # Check all bonds are single bonds
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
            return False, "Contains multiple bonds, not all bonds are single"

    # Check all carbons are sp3 hybridized (saturated)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon atom
            if atom.GetDegree() != 4:
                return False, "Carbon atom not saturated (not sp3 hybridized)"

    # Calculate molecular formula
    formula = rdMolDescriptors.CalcMolFormula(mol)
    # Parse formula to get counts of C and H
    from collections import Counter
    import re
    elements = re.findall('([A-Z][a-z]?)(\d*)', formula)
    element_counts = Counter()
    for (elem, count) in elements:
        count = int(count) if count else 1
        element_counts[elem] += count

    num_carbons = element_counts.get('C', 0)
    num_hydrogens = element_counts.get('H', 0)

    # Check general formula CnH2n+2
    if num_hydrogens != (2 * num_carbons + 2):
        return False, f"Does not fit general formula CnH2n+2 (found C{num_carbons}H{num_hydrogens})"

    return True, "Molecule is an acyclic saturated hydrocarbon fitting CnH2n+2"

__metadata__ = {   'chemical_class': {   'name': 'alkane',
                              'definition': 'An acyclic branched or unbranched hydrocarbon having the general formula CnH2n+2, consisting entirely of hydrogen atoms and saturated carbon atoms.'},
        'config': {   'llm_model_name': 'lbl/claude-sonnet'},
        'message': None,
        'attempt': 0,
        'success': True}