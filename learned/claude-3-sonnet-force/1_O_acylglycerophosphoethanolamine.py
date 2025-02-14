"""
Classifies: CHEBI:55493 1-O-acylglycerophosphoethanolamine
"""
"""
Classifies: CHEBI:85815 1-O-acylglycerophosphoethanolamine

A glycerophosphoethanolamine having an unspecified O-acyl substituent at the 1-position of the glycerol fragment.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_1_O_acylglycerophosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-O-acylglycerophosphoethanolamine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-O-acylglycerophosphoethanolamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerophosphoethanolamine backbone
    gpe_pattern = Chem.MolFromSmarts("[OX2][CHX4][CHX4][OX2][CX3](=[OX1])[OX2][CHX4][CHX4][OX2][P+](OCCN)(O)(O)")
    if not mol.HasSubstructMatch(gpe_pattern):
        return False, "No glycerophosphoethanolamine backbone found"
    
    # Look for acyl chain connected to the 1-position of the glycerol fragment
    acyl_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[OX2][CHX4][CHX4][OX2]")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern, maxMatches=1)
    
    if not acyl_matches:
        return False, "No acyl chain found"
    
    acyl_atom_idx = acyl_matches[0][-2]
    acyl_atom = mol.GetAtomWithIdx(acyl_atom_idx)
    
    if not any(bond.GetBeginAtomIdx() == acyl_atom_idx and bond.GetEndAtomIdx() in gpe_pattern.GetAtomsAsIndices() for bond in acyl_atom.GetBonds()):
        return False, "Acyl chain not connected to the 1-position of the glycerol fragment"

    # Count rotatable bonds to verify long acyl chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Acyl chain too short"

    # Check molecular weight - typically >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low"

    return True, "Contains glycerophosphoethanolamine backbone with acyl chain attached at the 1-position"