"""
Classifies: CHEBI:65323 monoterpenoid indole alkaloid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid indole alkaloid based on its SMILES string.
    A monoterpenoid indole alkaloid contains an indole moiety, a monoterpenoid unit and is an alkaloid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoterpenoid indole alkaloid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for a substituted Indole/Tryptophan-like substructure fused to a 6-member ring
    #   This pattern includes a substituted indole connected to a 6-membered ring (common in these alkaloids)
    #   The connection can be a direct bond or a spiro connection.
    #   The nitrogen can be substituted
    indole_fusion_pattern = Chem.MolFromSmarts("[c]1[c][nH][c]2[c]([c]1)[c]([#6])[c]([#6])[c]([#6])[#6]2")
    if not mol.HasSubstructMatch(indole_fusion_pattern):
       
        indole_fusion_pattern = Chem.MolFromSmarts("[c]1[c][N][c]2[c]([c]1)[c]([#6])[c]([#6])[c]([#6])[#6]2")
        if not mol.HasSubstructMatch(indole_fusion_pattern):
            
            indole_fusion_pattern = Chem.MolFromSmarts("[c]1[c][n][c]2[c]([c]1)[c]([#6])[c]([#6])[c]([#6])[#6]2") #nitrogen can be substituted
            if not mol.HasSubstructMatch(indole_fusion_pattern):

                indole_fusion_pattern = Chem.MolFromSmarts("[c]1[c][N+][c]2[c]([c]1)[c]([#6])[c]([#6])[c]([#6])[#6]2") #nitrogen can be positively charged
                if not mol.HasSubstructMatch(indole_fusion_pattern):

                    indole_fusion_pattern = Chem.MolFromSmarts("[c]1[c][nH+][c]2[c]([c]1)[c]([#6])[c]([#6])[c]([#6])[#6]2") #nitrogen can be protonated
                    if not mol.HasSubstructMatch(indole_fusion_pattern):
                        return False, "No appropriate substituted indole substructure fused to a ring found."


    # 2. Check for Complexity:
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings < 3: #increased threshold
        return False, "Too few rings for a monoterpenoid indole alkaloid."
        
    # Count carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 18: #increased threshold
        return False, "Too few carbons for a monoterpenoid indole alkaloid."

    # Count nitrogens
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 1:
         return False, "Must have at least one nitrogen (alkaloid nature)"

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300: # increased threshold
        return False, "Molecular weight is too low for a monoterpenoid indole alkaloid."

    return True, "Contains an appropriate indole substructure fused to a ring system, multiple rings, and sufficient complexity"