"""
Classifies: CHEBI:50699 oligosaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_oligosaccharide(smiles: str) -> tuple[bool, str]:
    """
    Determines if a molecule is an oligosaccharide based on its SMILES string.
    An oligosaccharide is a compound with monosaccharide units joined by glycosidic linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oligosaccharide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of multiple monosaccharide rings
    monosaccharide_rings = [r for r in mol.GetRingInfo().AtomRings() if len(r) in [5, 6]]
    if len(monosaccharide_rings) < 2:
        return False, "Fewer than two monosaccharide rings found"

    # Check for acetal oxygen atoms (indicative of glycosidic linkages)
    acetal_oxygens = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetHybridization() == Chem.HybridizationState.SP3 and atom.GetTotalDegree() == 2]

    # Check if acetal oxygens connect monosaccharide rings
    for o in acetal_oxygens:
        neighbors = [mol.GetAtomWithIdx(n).GetRingInfo().IsAtomInRingOfSize(5) or mol.GetAtomWithIdx(n).GetRingInfo().IsAtomInRingOfSize(6) for n in mol.GetAtomWithIdx(o).GetNeighbors()]
        if sum(neighbors) >= 2:
            break
    else:
        return False, "No glycosidic linkages found between monosaccharide rings"

    # Check molecular weight and composition
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for oligosaccharide"

    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if c_count < 10 or o_count < 5:
        return False, "Insufficient carbon or oxygen atoms for oligosaccharide"

    return True, "Contains multiple monosaccharide rings connected by glycosidic linkages"