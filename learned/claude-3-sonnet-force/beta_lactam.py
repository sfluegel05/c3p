"""
Classifies: CHEBI:35627 beta-lactam
"""
"""
Classifies: CHEBI:35458 beta-lactams
A lactam in which the amide bond is contained within a four-membered ring,
which includes the amide nitrogen and the carbonyl carbon.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
import re

def is_beta_lactam(smiles: str):
    """
    Determines if a molecule is a beta-lactam based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-lactam, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for 4-membered ring with N and C=O
    beta_lactam_pattern = Chem.MolFromSmarts("[NR2]1[CR2][CR2][CR2]1=O")
    if not mol.HasSubstructMatch(beta_lactam_pattern):
        return False, "No 4-membered ring with N and C=O found"

    # Check for common substituents and fused rings
    substituent_patterns = [
        Chem.MolFromSmarts("[OR1,NR1]"),  # OH, NH2, etc.
        Chem.MolFromSmarts("c1ccccc1"),   # Phenyl ring
        Chem.MolFromSmarts("C#N"),        # Nitrile
        Chem.MolFromSmarts("C=C"),        # Alkene
        Chem.MolFromSmarts("c1cccnc1"),   # Pyridine ring
        Chem.MolFromSmarts("c1ncccn1"),   # Pyrazine ring
        Chem.MolFromSmarts("c1nccnc1"),   # Pyrimidine ring
        Chem.MolFromSmarts("c1ncncc1"),   # Pyridazine ring
    ]
    has_common_substituent = any(mol.HasSubstructMatch(pattern) for pattern in substituent_patterns)

    # Check for ring strain and tautomerism
    ring_strain = rdMolDescriptors.CalcNumSpiroAtoms(mol) + rdMolDescriptors.CalcNumBridgeheadAtoms(mol)
    tautomer_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(mol, isomericSmiles=True), isomericSmiles=True))

    # Classify based on pattern, substituents, ring strain, and tautomerism
    if has_common_substituent and ring_strain > 0 and re.search(r"\[NH\]1\[CH\][CH\][CH\]1=O", tautomer_smiles):
        return True, "Contains strained 4-membered ring with N and C=O, and common beta-lactam substituents/fused rings"
    elif has_common_substituent and ring_strain > 0:
        return True, "Contains strained 4-membered ring with N and C=O, and common beta-lactam substituents/fused rings"
    elif ring_strain > 0 and re.search(r"\[NH\]1\[CH\][CH\][CH\]1=O", tautomer_smiles):
        return True, "Contains strained 4-membered ring with N and C=O (tautomeric form)"
    elif ring_strain > 0:
        return True, "Contains strained 4-membered ring with N and C=O"
    else:
        return False, "No characteristic features of beta-lactams found"