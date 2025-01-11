"""
Classifies: CHEBI:53156 polychlorobiphenyl
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polychlorobiphenyl(smiles: str):
    """
    Determines if a molecule is a polychlorobiphenyl (PCB) based on its SMILES string.
    PCBs have a biphenyl backbone with two or more chlorine substitutions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a PCB, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
    
    # Improved biphenyl detection: two benzene rings potentially connected in extended aromatic systems
    biphenyl_extended_pattern = Chem.MolFromSmarts("c1ccccc1-c2ccccc2")
    
    # Additional biphenyl derivatives, allowing for limited substitutions on benzene rings
    biphenyl_variants_pattern = Chem.MolFromSmarts("c1cc([cH])ccc1-c2c([cH])cc[H]c2")
    
    if not (mol.HasSubstructMatch(biphenyl_extended_pattern) or mol.HasSubstructMatch(biphenyl_variants_pattern)):
        return False, "No biphenyl backbone found or contains extensive substitutions inconsistent with PCBs"
    
    # Count the number of chlorine atoms
    cl_count = sum(atom.GetAtomicNum() == 17 for atom in mol.GetAtoms())
    
    if cl_count >= 2:
        # Identify and avoid additional complex functional groups not typical for standard PCBs
        functional_groups = rdMolDescriptors.CalcNumRings(mol) > 2  # Additional fused rings are unlikely in PCBs
        extraneous_oxygen_or_nitrogen = any(atom.GetAtomicNum() in [7, 8] for atom in mol.GetAtoms())
        
        if functional_groups or extraneous_oxygen_or_nitrogen:
            return False, "Structure is more complex than typical PCB with other atypical elements or functionality"

        return True, f"Contains biphenyl backbone with {cl_count} chlorine substitutions"
    else:
        return False, f"Contains biphenyl backbone but only {cl_count} chlorine substitutions, need at least 2"