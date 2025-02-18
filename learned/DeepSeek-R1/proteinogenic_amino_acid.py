"""
Classifies: CHEBI:83813 proteinogenic amino acid
"""
"""
Classifies: proteinogenic amino acid
"""
from rdkit import Chem
from rdkit.Chem import Mol
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import ChemicalFeatures
from rdkit import RDConfig
import os

def is_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a proteinogenic amino acid based on its SMILES string.
    These are alpha-amino acids with L-configuration (except glycine), including selenocysteine, pyrrolysine, and N-formylmethionine.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if proteinogenic amino acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES"
    
    # Check for alpha-amino acid structure: NH2-CH(R)-COOH
    # SMARTS pattern for N connected to alpha carbon connected to carboxyl
    # This allows for primary/secondary amines (proline) and different protonation states
    backbone_smarts = Chem.MolFromSmarts("[NX3]-[CX4H]([CX3](=O)[OX2H1])-*")
    if not mol.HasSubstructMatch(backbone_smarts):
        return False, "No alpha-amino acid backbone"
    
    # Check for at least one carboxylic acid group (redundant but safe)
    carboxyl_matches = mol.GetSubstructMatches(Chem.MolFromSmarts("[CX3](=O)[OX2H1]"))
    if not carboxyl_matches:
        return False, "No carboxylic acid group"
    
    # Find the alpha carbon(s) (connected to both N and COOH)
    alpha_carbons = set()
    for match in mol.GetSubstructMatches(backbone_smarts):
        alpha_carbon = match[1]
        alpha_carbons.add(alpha_carbon)
    
    # Check chirality for each alpha carbon (should be only one)
    for alpha_idx in alpha_carbons:
        alpha_atom = mol.GetAtomWithIdx(alpha_idx)
        chiral_tag = alpha_atom.GetChiralTag()
        
        # Glycine case: no chiral center, R is H
        if chiral_tag == Chem.ChiralType.CHI_UNSPECIFIED:
            # Check if R is H (only N and COOH attached)
            neighbors = alpha_atom.GetNeighbors()
            heavy_neighbors = [n for n in neighbors if n.GetAtomicNum() != 1]
            if len(heavy_neighbors) != 2:  # Expect N and COOH C
                return False, "Non-chiral alpha carbon with non-H substituents (not glycine)"
        else:
            # Check L-configuration (S in CIP if priorities are N > COOH > R > H)
            # Assign CIP labels
            Chem.AssignStereochemistry(mol)
            cip = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
            cip_dict = dict(cip)
            if alpha_idx not in cip_dict:
                return False, "Chiral center not assigned"
            
            cip_code = cip_dict[alpha_idx]
            if cip_code != 'S':
                return False, f"Alpha carbon has {cip_code} configuration, expected S (L)"
    
    return True, "Alpha-amino acid with correct backbone and L-configuration (if chiral)"