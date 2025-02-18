"""
Classifies: CHEBI:26377 pterocarpans
"""
"""
Classifies: CHEBI:51811 pterocarpan
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_pterocarpan(smiles: str):
    """
    Determines if a molecule is a pterocarpan based on its SMILES string.
    Pterocarpans are members of the class of benzofurochromene with a
    6a,11a-dihydro-6H-[1]benzofuro[3,2-c]chromene skeleton and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pterocarpan, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for pterocarpan scaffold pattern
    scaffold_pattern = Chem.MolFromSmarts("[C@]12Oc3ccccc3[C@@]1(O)COc1ccccc21")
    if not mol.HasSubstructMatch(scaffold_pattern):
        return False, "Does not contain the pterocarpan scaffold"

    # Check for fused rings
    fused_rings = AllChem.GetSymmSSSR(mol)
    if len(fused_rings) < 4:
        return False, "Does not have enough fused rings"

    # Check for required heteroatoms
    atom_counts = mol.GetAtomSMARTSPatternCount("O", 2, 4)
    if atom_counts <= 0:
        return False, "Does not contain the required oxygen atoms"

    # Check for common substitutions
    subs_pattern = Chem.MolFromSmarts("[OC]")
    subs_matches = mol.GetSubstructMatches(subs_pattern)
    if len(subs_matches) < 2:
        return False, "Does not have enough common substitutions"

    # Check for connectivity and size
    if not mol.GetRingInfo().AtomRings():
        return False, "Does not contain rings"
    if mol.GetNumHeavyAtoms() < 15:
        return False, "Too small to be a pterocarpan"

    return True, "Contains the pterocarpan scaffold with required fused rings and substitutions"