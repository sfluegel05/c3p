"""
Classifies: CHEBI:60027 polymer
"""
"""
Classifies: CHEBI:24702 polymer
A polymer is a mixture composed of macromolecules of different kinds, which may be differentiated by
composition, length, degree of branching, etc. This includes oligomers, cyclic polymers, and polymers
with different backbones or substituents.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_polymer(smiles: str):
    """
    Determines if a molecule is a polymer based on its SMILES string.

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

    # Check for common polymer backbones
    backbones = [
        Chem.MolFromSmarts("[C-]C(=O)O[C-]"),  # Esters
        Chem.MolFromSmarts("[C-]OC"),  # Ethers
        Chem.MolFromSmarts("[C-]C(=O)N[C-]"),  # Amides
        Chem.MolFromSmarts("[C-]C([C-])([C-])[C-]"),  # Alkanes
        Chem.MolFromSmarts("[C-]C([C-])=C([C-])[C-]"),  # Alkenes
        Chem.MolFromSmarts("[C-]#C[C-]"),  # Alkynes
        Chem.MolFromSmarts("[C-]C([C-])([C-])[Si-]"),  # Siloxanes
        Chem.MolFromSmarts("[N-][C-]=[N-]"),  # Azides
    ]
    if not any(mol.HasSubstructMatch(backbone) for backbone in backbones):
        return False, "No common polymer backbone found"

    # Check for repeating units
    sssr = Chem.GetSymmSSSR(mol)
    if not sssr:
        return False, "No repeating units found"

    # Check for molecular weight (flexible threshold for oligomers)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for polymer"

    # Check for degree of branching
    branching = rdMolDescriptors.CalcDegreeOfBranching(mol)
    if branching < 0.2:
        return True, "Linear polymer"
    elif branching < 0.8:
        return True, "Branched polymer"
    else:
        return True, "Highly branched polymer"

    # If all criteria are met, classify as polymer
    return True, "Molecule exhibits characteristics of a polymer"