"""
Classifies: CHEBI:26607 saturated fatty acid
"""
"""
Classifies: CHEBI:36713 saturated fatty acid
Any fatty acid containing no carbon to carbon multiple bonds. Known to produce adverse biological effects when ingested to excess.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a saturated fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a saturated fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles, removeHs=False)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxyl group (-C(=O)O)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[O;H,-]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxyl group found"

    # Check for absence of multiple bonds
    if any(bond.GetBondType() == Chem.BondType.DOUBLE for bond in mol.GetBonds()):
        return False, "Contains carbon-carbon double bonds"
    if any(bond.GetBondType() == Chem.BondType.TRIPLE for bond in mol.GetBonds()):
        return False, "Contains carbon-carbon triple bonds"

    # Check for aliphatic carbon chain (allow branching)
    aliphatic_pattern = Chem.MolFromSmarts("C~C~C~C")  # Recursive pattern for variable-length chain
    aliphatic_matches = mol.GetSubstructMatches(aliphatic_pattern)
    if not aliphatic_matches:
        return False, "No aliphatic carbon chain found"

    # Check for branching patterns
    branched_pattern = Chem.MolFromSmarts("[C;D4]")  # Quaternary carbon indicates branching
    branched_matches = mol.GetSubstructMatches(branched_pattern)

    # Check molecular weight (at least 120 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 120:
        return False, "Molecular weight too low for fatty acid"

    # Check for additional functional groups (optional)
    alcohol_pattern = Chem.MolFromSmarts("[OX1H]")
    alcohol_matches = mol.GetSubstructMatches(alcohol_pattern)
    methyl_pattern = Chem.MolFromSmarts("C[C;D4]")
    methyl_matches = mol.GetSubstructMatches(methyl_pattern)

    if branched_matches:
        if alcohol_matches:
            return True, "Contains a carboxyl group, a saturated branched aliphatic carbon chain, and an alcohol group"
        elif methyl_matches:
            return True, "Contains a carboxyl group, a saturated branched aliphatic carbon chain, and methyl groups"
        else:
            return True, "Contains a carboxyl group and a saturated branched aliphatic carbon chain"
    else:
        if alcohol_matches:
            return True, "Contains a carboxyl group, a saturated linear aliphatic carbon chain, and an alcohol group"
        elif methyl_matches:
            return True, "Contains a carboxyl group, a saturated linear aliphatic carbon chain, and methyl groups"
        else:
            return True, "Contains a carboxyl group and a saturated linear aliphatic carbon chain"