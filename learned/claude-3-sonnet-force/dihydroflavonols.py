"""
Classifies: CHEBI:48039 dihydroflavonols
"""
"""
Classifies: CHEBI:33716 dihydroflavonol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_dihydroflavonols(smiles: str):
    """
    Determines if a molecule is a dihydroflavonol based on its SMILES string.
    A dihydroflavonol is a hydroxyflavanone with a hydroxy group at position 3 of the heterocyclic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dihydroflavonol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for the dihydroflavonol scaffold
    dihydroflavonol_pattern = Chem.MolFromSmarts("[O;R]1[C@H]([C@@H](C(=O)C2=C1C=CC=C2)O)[c;R]2ccccc2")
    if not mol.HasSubstructMatch(dihydroflavonol_pattern):
        return False, "Molecule does not contain the dihydroflavonol scaffold"
    
    # Check for hydroxy group at position 3 of the heterocyclic ring
    hydroxy_at_pos_3 = any(atom.GetSymbol() == 'O' and atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3
                           for atom in mol.GetAtomWithIdx(2).GetNeighbors())
    if not hydroxy_at_pos_3:
        return False, "Molecule does not have a hydroxy group at position 3 of the heterocyclic ring"
    
    # Count the number of hydroxy groups
    hydroxyl_groups = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'O' and atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3)
    if hydroxyl_groups < 3:
        return False, "Molecule does not have enough hydroxy groups"
    
    # Check molecular weight - dihydroflavonols typically 200-600 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 600:
        return False, "Molecular weight outside typical range for dihydroflavonols"
    
    return True, "Molecule meets the criteria for a dihydroflavonol"