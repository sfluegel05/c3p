"""
Classifies: CHEBI:72588 semisynthetic derivative
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_semisynthetic_derivative(smiles: str):
    """
    Determines if a molecule is a semisynthetic derivative based on its SMILES string.
    This version focuses on identifying a complex scaffold and synthetic modifications

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

    if num_heavy_atoms < 15:
        return False, "Too few heavy atoms to be a semisynthetic derivative."

    if num_rings < 1 and num_chiral_centers < 1:
        return False, "Too simple, lacks ring structures and chiral centers"

    # Check for complex scaffold
    # (This is a very basic check, ideally should use a database of SMARTS patterns)
    if num_rings < 2 and num_chiral_centers < 3:
        return False, "Molecule lacks a complex scaffold, unlikely a semi-synthetic."
    
    # Common synthetic modifications (esters, amides, ethers, halogens)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    num_esters = len(mol.GetSubstructMatches(ester_pattern))
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")
    num_amides = len(mol.GetSubstructMatches(amide_pattern))
    ether_pattern = Chem.MolFromSmarts("[OX2][CX4]")
    num_ethers = len(mol.GetSubstructMatches(ether_pattern))
    num_halogens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() in [9, 17, 35, 53])
    
    #Check for combination of natural core with synthetic modifications. 
    if (num_esters + num_amides + num_ethers) > 0 or num_halogens > 0:
            return True, "Contains structural features suggestive of a semisynthetic derivative"
    else:
        return False, "No common synthetic modifications found"