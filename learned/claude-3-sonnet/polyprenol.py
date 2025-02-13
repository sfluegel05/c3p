"""
Classifies: CHEBI:26199 polyprenol
"""
"""
Classifies: CHEBI:35641 polyprenol

A polyprenol is defined as: Any member of the class of prenols possessing the general formula
H-[CH2C(Me)=CHCH2]nOH in which the carbon skeleton is composed of more than one isoprene units.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyprenol(smiles: str):
    """
    Determines if a molecule is a polyprenol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyprenol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for isoprene unit pattern [CH2]=[CH][C@H]([CH3])[CH2]=[CH2]
    isoprene_pattern = Chem.MolFromSmarts("[CH2]=[CH][C@H]([CH3])[CH2]=[CH2]")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) < 2:
        return False, "Less than two isoprene units found"

    # Check for terminal -OH group
    has_terminal_oh = any(atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 1 for atom in mol.GetAtoms())
    if not has_terminal_oh:
        return False, "No terminal -OH group found"

    # Check for linear or branched carbon chain
    chain_pattern = Chem.MolFromSmarts("[CH2,CH3][CH2,CH3]~[CH2,CH3]~[CH2,CH3]~[CH2,CH3]")
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    if not chain_matches:
        return False, "No linear or branched carbon chain found"

    # Check molecular weight range (typical polyprenols are between 200-1000 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 1000:
        return False, "Molecular weight outside typical range for polyprenols"

    return True, "Molecule contains multiple isoprene units, a terminal -OH group, and a linear or branched carbon chain"