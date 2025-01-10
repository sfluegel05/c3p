"""
Classifies: CHEBI:29017 1-acyl-sn-glycero-3-phosphoethanolamine
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_1_acyl_sn_glycero_3_phosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine based on its SMILES string.
    This class is characterized as a 1-O-acylglycerophosphoethanolamine with (R)-configuration at sn-2.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Recognize the glycerol backbone with sn-2 chiral center
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    # Look for (R) configuration; this is complex, need structural context
    sn2_chiral_center = None
    for center in chiral_centers:
        atom_idx, chirality = center
        atom = mol.GetAtomWithIdx(atom_idx)
        # Assume that C bonded to phosphate group and secondary alcohol as sn-2
        neighbors = [neigh.GetSymbol() for neigh in atom.GetNeighbors()]
        if 'O' in neighbors and 'P' in neighbors:
            sn2_chiral_center = (atom_idx, chirality)
            break

    if sn2_chiral_center is None or sn2_chiral_center[1] != 'R':
        return False, "No glycerol backbone with (R)-configuration at sn-2 position found"

    # Check for the glycerol-like structure with ester at sn-1
    glycerol_sn1_ester_pattern = Chem.MolFromSmarts("O[C@H](CO)C(=O)C")
    if not mol.HasSubstructMatch(glycerol_sn1_ester_pattern):
        return False, "No ester linkage at sn-1 position detected"

    # Confirm the attached phosphoethanolamine group
    phosphate_ethanolamine_pattern = Chem.MolFromSmarts("P(=O)([O-])[O]OCCN")
    if not mol.HasSubstructMatch(phosphate_ethanolamine_pattern):
        return False, "No phosphoethanolamine group found"
    
    return True, "Molecule matches 1-acyl-sn-glycero-3-phosphoethanolamine structural criteria"