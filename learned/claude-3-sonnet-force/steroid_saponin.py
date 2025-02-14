"""
Classifies: CHEBI:61655 steroid saponin
"""
"""
Classifies: CHEBI:28304 steroid saponin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_steroid_saponin(smiles: str):
    """
    Determines if a molecule is a steroid saponin based on its SMILES string.
    A steroid saponin is a saponin derived from a hydroxysteroid backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a steroid saponin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for steroid backbone patterns
    steroid_patterns = [
        Chem.MolFromSmarts("[C@H]1CC[C@]2(C)[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@@]21C"),  # Androstane
        Chem.MolFromSmarts("[C@H]1CC[C@]2(C)[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@@]12C"),  # Pregnane
        Chem.MolFromSmarts("[C@H]1CC[C@]2(C)[C@@H]3C[C@H](O)C=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@@]12C"),  # Cholane
        Chem.MolFromSmarts("[C@H]1C[C@@]2([C@@H](C[C@@]3([C@@H](C[C@@]4([C@@H](C[C@@H](O[C@@H]5O[C@@H]([C@@H](O)[C@H](O)[C@H]5O)CO)C4)C)[H])C=C3)C)(C2)C)[H]")  # Spirostane
        # Add more steroid backbone patterns as needed
    ]
    has_steroid_backbone = any(mol.HasSubstructMatch(pattern) for pattern in steroid_patterns)
    if not has_steroid_backbone:
        # Check for other potential steroid backbones or rearranged skeletons
        steroid_features = [
            Chem.MolFromSmarts("[C@H]1CC[C@]2(C)[C@@H]3[C@H](C[C@H]4[C@@]5(CC[C@@H](C6=CC(=O)CC6)[C@@H](C5)CC4)C)[C@H]3CC[C@@]12C"),  # Cardanolide
            Chem.MolFromSmarts("[C@H]1CC[C@]2(C)[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3C[C@@]12C")  # Estrane
            # Add more patterns for other steroid backbones
        ]
        has_steroid_features = any(mol.HasSubstructMatch(pattern) for pattern in steroid_features)
        if not has_steroid_features:
            return False, "No steroid backbone found"
    
    # Look for glycosidation patterns
    sugar_patterns = [
        Chem.MolFromSmarts("[OX2][CX4][CX4][OX2]"),  # Glycosidic linkage
        Chem.MolFromSmarts("[OX2][CX3](=[OX1])[CX4][CX4][OX2]"),  # Glycoside group
        Chem.MolFromSmarts("C1CCCCC1O"),  # Hexose ring
        Chem.MolFromSmarts("C1CCCC1O")  # Pentose ring
        # Add more patterns for sugar moieties or glycosidation
    ]
    has_glycosidation = any(mol.HasSubstructMatch(pattern) for pattern in sugar_patterns)
    if not has_glycosidation:
        return False, "No evidence of glycosidation found"
    
    # Check for additional structural features (optional)
    has_other_rings = bool(Chem.GetSymmSSSR(mol))
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    # Apply additional rules or heuristics as needed (optional)
    if mol_wt < 400:
        return False, "Molecular weight too low for typical steroid saponin"
    if n_rotatable < 3:
        return False, "Too few rotatable bonds for typical steroid saponin"
    
    return True, "Contains steroid backbone and evidence of glycosidation, consistent with a steroid saponin structure"