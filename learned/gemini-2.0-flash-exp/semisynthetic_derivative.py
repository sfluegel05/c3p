"""
Classifies: CHEBI:72588 semisynthetic derivative
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_semisynthetic_derivative(smiles: str):
    """
    Determines if a molecule is a semisynthetic derivative based on its SMILES string.
    This version focuses on identifying a natural product-like core and synthetic modifications.
    It uses SMARTS patterns for common natural product cores and synthetic modifications.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely a semisynthetic derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    num_heavy_atoms = mol.GetNumHeavyAtoms()
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    num_chiral_centers = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
    
    if num_heavy_atoms < 10: #Added min heavy atoms
        return False, "Too few heavy atoms to be a semisynthetic derivative."
    
    if num_rings <= 1 and num_chiral_centers < 2:
        #Initial check to filter linear/simple structures
        return False, "Too simple, lacks ring structures and chiral centers"

    # 1. Natural product core detection (more specific than before)
    natural_product_cores = [
        Chem.MolFromSmarts("[C]1[C]([C])([C])[C][C]1[C]"), #terpene core
        Chem.MolFromSmarts("[C]12[C]([C])([C])[C][C]1[C]3[C][C][C][C]23"), #fused terpene
        Chem.MolFromSmarts("[C]12[C]([C])([C])[C][C]1[C][C]3[C][C][C][C]23"), #fused terpene
        Chem.MolFromSmarts("[C]12[C]([C])([C])[C][C]1[C][C][C][C]2"), #fused terpene with a 4 membered ring
        Chem.MolFromSmarts("[C]1[C]([C])([C])[C][C]1[C]2[C][C][C]2"), #fused terpene
        Chem.MolFromSmarts("[C]1[C]2[C]([C]([C])([C])[C][C]2)[C][C]1[C]"), #fused terpene another
        Chem.MolFromSmarts("[C]1[C][C]2[C]([C])([C])[C][C]2[C][C]1"), #fused terpene more
        Chem.MolFromSmarts("[C]1[C][C]2[C]([C])([C])[C][C]2[C][C][C]1"), #fused terpene more
        Chem.MolFromSmarts("[C]12[C]([C])([C])[C][C]1[C][C]3[C][C][C]32"), #fused terpene with 2 5 membered rings
        Chem.MolFromSmarts("[C]12[C]([C])([C])[C][C]1[C]3[C][C][C]23"), #fused terpene with fused 5 and 6 membered rings
        Chem.MolFromSmarts("[C]12[C]([C])([C])[C][C]1[C][C][C][C][C]2"), #fused terpene with a six membered ring
        Chem.MolFromSmarts("[C]12[C]([C])([C])[C][C]1[C][C]3[C][C][C][C]32"), #fused terpene another 6 membered ring
        Chem.MolFromSmarts("[C]12[C]([C])([C])[C][C]1[C]3[C][C][C][C]32"), #fused terpene another 6 membered ring
        Chem.MolFromSmarts("C[C@]12[C@@H]([C@@H]([C@H](CC[C@@H]1CC[C@@H]3[C@@H]2CC=C4[C@@H]3[C@@H](C=C[C@@H]5[C@@H]4CC[C@H](C5)O)C)C)C)"), #Steroid Core
        Chem.MolFromSmarts("[C]1[C]([C])([C])[C][C]1[C]2[C][C][C][C]2[C]3[C][C][C]3"), #fused rings like in steroids
        Chem.MolFromSmarts("[C]1[C]([C])([C])[C][C]1[C]2[C][C][C][C]2[C]3[C][C][C][C]3"), #fused rings like in steroids more
        Chem.MolFromSmarts("[C]12[C]([C])([C])[C][C]1[C][C]3[C][C][C][C]32"), #fused rings like in steroids
        Chem.MolFromSmarts("[C]1[C]([C])([C])[C][C]1[C][C]2[C][C][C]2[C]3[C][C][C][C]3"), #fused rings like in steroids
        Chem.MolFromSmarts("[C]12[C]([C])([C])[C][C]1[C][C]3[C][C][C]32"), #fused rings like in steroids another
        Chem.MolFromSmarts("C[C@]1([C@@H]([C@@H]([C@H]([C@@H](O1)O[C@@H]2[C@@H]([C@H]([C@@H]([C@@H](O2)C)O)OC)C)OC)C)"),  #Macrolide Core
        Chem.MolFromSmarts("[C]1[C]([C])([C])[C][C]1[C]2[C][C]3[C][C][C][C]32") # Core similar to ergot alkaloids
    ]
    has_natural_core = any(mol.HasSubstructMatch(core) for core in natural_product_cores)


    if not has_natural_core:
        return False, "Molecule does not contain a recognizable natural product-like core."

    # 2. Synthetic modification detection (more specific)
    synthetic_modifications = [
        Chem.MolFromSmarts("[CX4]([Cl,Br,I,F])"),   # Halogenated carbons
        Chem.MolFromSmarts("[CX3](=[OX1])[OX2][CX4](C)C"), #Acetate
        Chem.MolFromSmarts("[CX3](=[OX1])[OX2][CX4][CX4]"), # Generic ester
        Chem.MolFromSmarts("[NX3][CX3](=[OX1])[CX4]"), # Amide next to a carbon
        Chem.MolFromSmarts("[NX3][CX3](=[OX1])[C;!H0]"), # amide, not formamide
        Chem.MolFromSmarts("[OX2][CX4][OX2]"), # Ether
        Chem.MolFromSmarts("[CX3](=[OX1])[OX2][CX3](=[OX1])"), # Anhydride
        Chem.MolFromSmarts("[C]#[C][C]"), # Alkyne
        Chem.MolFromSmarts("[C]=[C]=[C]"),#Allene
        Chem.MolFromSmarts("[C]1[C]=[C][C]=[C][C]1"), #1,3,5-substituted 6 membered ring
        Chem.MolFromSmarts("[C]1[C]=[C][C][C]=[C]1"), #1,4-substituted 6 membered ring
        Chem.MolFromSmarts("[C]1[C][C]=[C][C][C]1"), # 1,2-substituted 6 membered ring
        Chem.MolFromSmarts("[C]1[C]([C])=[C][C]([C])=[C]1"), # substituted aromatic ring
        Chem.MolFromSmarts("[C]1[C]([C])=[C][C]([C])=[C][C]1"),# substituted aromatic ring more
        Chem.MolFromSmarts("c1ccccc1[C;!H0]"), # substituted benzene
        Chem.MolFromSmarts("[C]1([C])=[C][C]=[C][C]([C])=[C]1"), # Substituted aromatic ring
        Chem.MolFromSmarts("[C]1=[C][C]([C])=[C][C]=[C]1"), # Substituted aromatic ring more
        Chem.MolFromSmarts("c1c([C])cc([C])c([C])c1"),#trisubstituted benzene
        Chem.MolFromSmarts("[C;!H0][N]([C;!H0])=[C][N]"), # substituted amide or urea
        Chem.MolFromSmarts("[C;!H0][N]([C;!H0])=[C][O]"), # substituted carbamate
    ]
   
    modification_count = sum(len(mol.GetSubstructMatches(mod)) for mod in synthetic_modifications)


    if modification_count >= 2:
        return True, "Contains structural features suggestive of a semisynthetic derivative"
    else:
         return False, "Contains a natural product core, but with too few synthetic modifications"