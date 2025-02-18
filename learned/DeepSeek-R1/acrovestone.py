"""
Classifies: CHEBI:2440 acrovestone
"""
"""
Classifies: acrovestone (a polyphenol from Acronychia pedunculata with antioxidant/antityrosinase activity)
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_acrovestone(smiles: str):
    """
    Determines if a molecule is acrovestone based on its SMILES string.
    Acrovestone is characterized as an isoflavone derivative with glycosylation and multiple phenolic groups.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if matches acrovestone characteristics, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for isoflavone core (benzopyran-4-one with phenyl substituent)
    isoflavone_pattern = Chem.MolFromSmarts("[#6]1=[#8]-[#6]2:[#6]:[#6]:[#6](:[#6]:[#6]:1)-[#6](=[#8])-[#6]2-[#6]~[#6]~[#6]~[#6]~[#6]~[#6]")
    if not mol.HasSubstructMatch(isoflavone_pattern):
        return False, "No isoflavone core detected"

    # Check for glycosylation (O-linked to carbohydrate-like structure)
    glycoside_pattern = Chem.MolFromSmarts("[OX2]-[#6;r5,r6]")
    if not mol.HasSubstructMatch(glycoside_pattern):
        return False, "No glycosylation detected"

    # Verify polyphenol characteristics (multiple phenolic -OH groups)
    phenolic_oh = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 1:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetHybridization() == Chem.HybridizationType.SP2:
                    phenolic_oh += 1
                    break
    if phenolic_oh < 2:
        return False, f"Insufficient phenolic OH groups ({phenolic_oh})"

    # Check molecular weight range based on examples (300-1000 Da)
    mol_wt = Descriptors.ExactMolWt(mol)
    if not (300 < mol_wt < 1000):
        return False, f"Molecular weight {mol_wt:.1f} outside typical range"

    return True, "Isoflavone core with glycosylation and phenolic groups"